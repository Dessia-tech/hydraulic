#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 14:24:59 2018

@author: steven
"""

import math
import numpy as npy
import matplotlib.pyplot as plt
import dessia_common as dc
import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.primitives3d as p3d
from copy import copy

from typing import List

# Definition of equivalent L/D values
# Lists in increasing order
ben_angle = [0, 45, 90, 180, 360]
ben_LD = [0, 16, 30, 60, 120]
# Bend radius ?
mbd_angle = [0, 45, 90]
mbd_LD = [0, 15, 58]
enl_rap = [1/4, 1/2, 3/4, 1]
enl_LD = [30, 20, 6.5, 0]
ctr_rap = [1, 4/3, 2, 4]
ctr_LD = [0, 6.5, 11, 15]

class StraightPipe(dc.DessiaObject):
    """
    Abstract Class
    Straight pipes linking 2 points
    """
    def __init__(self, point1, point2, diameter:float,
                 heat_exchange:bool=True, name:str=''):
        self.points = [point1, point2]
        self.active_points = self.points
        self.heat_exchange = heat_exchange
        self.radius = diameter/2
        self.surf = math.pi*self.radius**2
        self.length = point1.point_distance(point2)
        self.fQ = 16*self.length/(math.pi*self.radius**4)
        self.n_equations = 2
        self.name = name

    def __str__(self):
        # return "{} from {} to {}".format(self.__class__.__name__, self.points[0], self.points[1])
        return 'Straight pipe'

    def system_matrix(self, constant):
        system_matrix = npy.array([[-constant, self.fQ, constant, 0],
                                   [0, 1, 0, 1]])
        return system_matrix

    def gmsh_repr1d(self, j, points_index):
        line = "Line({}) = [{},{}];\n".format(j, points_index[self.points[0]],
                                              points_index[self.points[1]])
        return (line, [j])


class StraightPipe2D(StraightPipe):
    """
    Straight pipes linking 2 2D points
    """
    def __init__(self, point1:vm.Point2D, point2:vm.Point2D, diameter:float,
                 heat_exchange:bool=True, name:str=''):
        StraightPipe.__init__(self, point1, point2, diameter, heat_exchange,
                              name=name)

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        ax = vme.LineSegment2D(*self.points).plot(ax)
        return ax

class StraightPipe3D(StraightPipe):
    """
    Straight pipes linking 2 3D points
    """
    def __init__(self, point1:vm.Point3D, point2:vm.Point3D, diameter:float,
                 heat_exchange:bool=True, name:str=''):
        StraightPipe.__init__(self, point1, point2, diameter, heat_exchange,
                              name=name)

    def plot2d(self, x3D, y3D, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        vme.LineSegment3D(*self.points).plot2d(x3D, y3D, ax)
        return ax

    def plot(self, ax=None):
        ax = vme.LineSegment3D(*self.points).plot(ax=ax)
        return ax

    def volmdlr_primitives(self):
        axis = self.points[1] - self.points[0]
        axis.normalize()
        return [p3d.HollowCylinder(0.5*(self.points[0] + self.points[1]), axis,
                                  self.radius, self.radius + 0.001,
                                  self.length,
                                  name=self.name)]


class SingularPipe(dc.DessiaObject):
    """
    Other type of pipes linking 2 points
    """
    def __init__(self, p1, p2, form, heat_exchange=True, name=''):
        self.points = [p1, p2]
        self.active_points = self.points
        self.type = form
        self.heat_exchange = heat_exchange
        self.name = name

    def __str__(self):
        return "{}-{}-{}".format(self.points[0], self.type, self.points[1])

    def system_matrix(self, constant):
        system_matrix = npy.array([[-constant, self.fQ, constant, 0],
                                   [0, 1, 0, 1]])
        return system_matrix

class Bend(SingularPipe):
    """
    Bend between two points, defined by third, interior point
    """
    def __init__(self, start_point, interior_point, end_point,
                 arc, diameter,
                 heat_exchange=True,
                 name=''):
        SingularPipe.__init__(self, start_point, end_point, 'ben', heat_exchange, name)
        self.start_point = start_point
        self.interior_point = interior_point
        self.end_point = end_point
        self.radius = 0.5 * diameter
        self.arc = arc
        self.turn_radius = self.arc.radius
        self.turn_angle = self.arc.angle
        self.section = math.pi*self.radius**2
        self.length = self.turn_radius*self.turn_angle
        length_d = get_equivalent(abs(self.turn_angle*180/math.pi), ben_LD, ben_angle)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)
        self.n_equations = 2

    def gmsh_repr1d(self, j, points_index):
        # Returns the 1D representation of the pipe for gmsh
        n = max(points_index.values())+1
        circle = "Point({0}) = [{1[0]},{1[1]},0.,1.];\n".format(n+1, self.arc.center)
        points_index["center of {}".format(str(self))] = n+1
        circle += "Point({0}) = [{1[0]},{1[1]},0.,1.];\n".format(n+2, self.interior_point)
        points_index["intermediate point of {}".format(str(self))] = n+2
        circle += "Circle({}) = [{},{},{}];\n".format(j, points_index[self.points[0]], n+1, n+2)
        circle += "Circle({}) = [{},{},{}];\n".format(j+1, n+2, n+1, points_index[self.points[1]])
        phys = [j, j+1]
        return (circle, phys)


class Bend2D(Bend):
    """
    Bend between two points, defined by third, interior point
    """
    def __init__(self, start_point, interior_point, end_point,
                 diameter,
                 heat_exchange=True,
                 name=''):
        arc = vme.Arc2D(start_point, interior_point, end_point)
        Bend.__init__(self, start_point, interior_point, end_point,
                      arc, diameter, heat_exchange, name)

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        self.arc.plot(ax)

class Bend3D(Bend):
    """
    Bend between two points, defined by third, interior point
    """
    def __init__(self, start_point, interior_point, end_point,
                 diameter, heat_exchange=True, name=''):
        arc = vme.Arc3D(start_point, interior_point, end_point)
        Bend.__init__(self, start_point, interior_point, end_point,
                      arc, diameter, heat_exchange, name)

    def plot2d(self, x3D, y3D, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        self.arc.plot2d(x3D, y3D, ax=ax)
        return ax
        
    def plot(self, ax=None):
        self.arc.plot(ax=ax)
        return ax

    def volmdlr_primitives(self):
        # normal_section = (self.arc.start - self.arc.center).cross(self.arc.normal)
        section = vmw.Circle2D(vm.O2D, self.radius+0.001)
        return [p3d.Sweep(section, vmw.Wire3D([self.arc]))]


class MitterBend(SingularPipe):
    """
    Mitter bend
    (For one point, specify 2 points with same coordinates)
    """
    def __init__(self, p1, p2, diameter, angle, heat_exchange=True, name=''):
        SingularPipe.__init__(self, p1, p2, 'mbd', heat_exchange, name)
        self.radius = diameter/2
        self.surf = math.pi*self.radius**2
        self.length = p1.PointDistance(p2)
        length_d = get_equivalent(angle*180/math.pi, mbd_LD, mbd_angle)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)

    def gmsh_repr1D(self, j, points_index):
        # Returns the 1D representation of the pipe for gmsh
        line = "Line({}) = [{},{}];\n".format(j, points_index[self.points[0]],
                                              points_index[self.points[1]])
        return (line, [j])

class Enlargement(SingularPipe):
    """
    Enlargement from diameter 1 at point 1 to diameter 2 at point 2
    for gradual enlargement, specify the angle
    """
    def __init__(self, p1, p2, d1, d2, angle=None, heat_exchange=True, name=''):
        SingularPipe.__init__(self, p1, p2, 'enl', heat_exchange, name)
        self.radius = d1/2
        self.surf = math.pi*(d1/2)**2
        self.length = p1.PointDistance(p2)
        length_d = get_equivalent(d1/d2, enl_LD, enl_rap)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)

    def gmsh_repr1d(self, j, points_index):
        # Returns the 1D representation of the pipe for gmsh
        line = "Line({}) = [{},{}];\n".format(j, points_index[self.points[0]],
                                              points_index[self.points[1]])
        return (line, [j])

class Contraction(SingularPipe):
    """
    Contraction from diameter 1 at point 1 to diameter 2 at point 2
    for gradual ontraction, specify the angle
    """
    def __init__(self, p1, p2, d1, d2, angle=None, heat_exchange=True, name=''):
        SingularPipe.__init__(self, p1, p2, 'ctr', heat_exchange, name)
        self.radius = d1/2
        self.surf = math.pi*(d1/2)**2
        self.length = p1.PointDistance(p2)
        length_d = get_equivalent(d1/d2, ctr_LD, ctr_rap)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)

    def gmsh_repr1d(self, j, points_index):
        # Returns the 1D representation of the pipe for gmsh
        line = "Line({}) = [{},{}];\n".format(j, points_index[self.points[0]],
                                              points_index[self.points[1]])
        return (line, [j])

#class UserDefined(SingularPipe):
#    """
#    Singular Pipe defined by user with a known parameter and different
#    L/D equivalent ratios for different values of the parameter.
#    """
#    def __init__(self, p1, p2,
#                 diameter,
#                 param_values,
#                 LD_values,
#                 *params,
#                 heat_exchange=True,
#                 name=''):
#        SingularPipe.__init__(self, p1, p2, 'usr', heat_exchange, name)
#        self.radius = diameter/2
#        self.surf = math.pi*self.radius**2
#        self.length = p1.PointDistance(p2)
#        length_d = GetEquivalent(param_values, LD_values, *params)
#        self.fQ = 16*2*length_d/(math.pi*self.radius**3)
#
#    def gmsh_repr1d(self, j, points_index):
#        # Returns the 1D representation of the pipe for gmsh
#        line = "Line({}) = [{},{}];\n".format(j, points_index[self.points[0]],
#                                              points_index[self.points[1]])
#        return (line, [j])

class UserDefined(SingularPipe):
    """
    Singular Pipe defined by user with a known parameter and different
    L/D equivalent ratios for different values of the parameter.
    """
    def __init__(self, p1, p2,
                 fQ,
                 heat_exchange=True,
                 name=''):
        SingularPipe.__init__(self, p1, p2, 'usr', heat_exchange, name)
        self.fQ = fQ
        self.n_equations = 2
        
    def system_matrix(self, constant):
        system_matrix = npy.array([[-constant, self.fQ, constant, 0],
                                   [0, 1, 0, 1]])
        return system_matrix
        
    def volmdlr_primitives(self):
        axis = self.points[1] - self.points[0]
        l = axis.Norm()
        axis.Normalize()
        
        return [p3d.HollowCylinder(0.5*(self.points[0] + self.points[1]), axis,
                                  0.001, 0.002,
                                  l,
                                  name=self.name)]

    def gmsh_repr1d(self, j, points_index):
        # Returns the 1D representation of the pipe for gmsh
        line = "Line({}) = [{},{}];\n".format(j, points_index[self.points[0]],
                                              points_index[self.points[1]])
        return (line, [j])
    
    def plot(self, x3D, y3D, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        vm.LineSegment3D(*self.points).plot2D(x3D, y3D, ax)

class JunctionPipe(dc.DessiaObject):
    """
    Add pressure drop values
    junction linking 1 pipe to 2+ pipes
    """
    def __init__(self, central_point:vm.Point3D, other_points:List[vm.Point3D],
                 diameter:float, name:str=''):
        self.points = [central_point] + other_points
        self.active_points = other_points
        self.central_point = central_point
        self.heat_exchange = False
        self.radius = diameter/2
        self.surf = math.pi*self.radius**2
        self.lengths = [point.point_distance(central_point) for point in other_points]
        self.n_equations = 1 + len(self.active_points)
        self.name = name

        # To be replaced by values for junction pipes
        length_d = [length/diameter for length in self.lengths]
        self.fQs = [16*2*ld/(math.pi*self.radius**3) for ld in length_d]


    def __str__(self):
        # return "{}-jun-{}".format(self.central_point, len(self.active_points))
        return 'Junction'

    def system_matrix(self, constant):
        """
        Build local matrices
        """
        matrix_generator = [[0, 0] + [0, 1]*len(self.active_points)]
        for i in range(len(self.active_points)):
            i_real = i+1
            fq_value = self.fQs[i]
            line_generator = [0]*len(self.points)*2
            line_generator[0] = -constant
            line_generator[2*i_real] = constant
            line_generator[2*i_real + 1] = -fq_value
            matrix_generator.append(line_generator)

        system_matrix = npy.array(matrix_generator)
        return system_matrix

    def plot2d(self, x3D=vm.X3D, y3D=vm.Y3D, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        for point in self.active_points:
            vme.LineSegment2D(point, self.central_point).plot(ax)
        return ax

    def plot(self, ax=None):
        
        ax = vme.LineSegment3D(self.active_points[0], self.central_point).plot(ax=ax)
        
        for point in self.active_points[1:]:
            vme.LineSegment3D(point, self.central_point).plot(ax)
        return ax

    def gmsh_repr1d(self, j, points_index):
        # Returns the 1D representation of the pipe for gmsh
        line = "Line({}) = [{},{}];\n".format(j, points_index[self.active_points[0]],
                                              points_index[self.active_points[1]])
        return(line, [j])
    
    def volmdlr_primitives(self):
        primitives = []
        for point in self.active_points:
            axis = point - self.central_point
            axis.normalize()
            length = point.point_distance(self.central_point)
            primitives.append(p3d.HollowCylinder(0.5*(self.central_point + point), axis,
                                                 self.radius, self.radius + 0.001,
                                                 length,
                                                 name=self.name))
        return primitives

VM_EQUIVALENCES = {vme.LineSegment2D: (StraightPipe2D, (('start', None), ('end', None))),
                   vme.Arc2D: (Bend2D, (('start', None), ('interior', None), ('end', None))),
                   vme.LineSegment3D: (StraightPipe3D, (('start', None), ('end', None))),
                   vme.Arc3D: (Bend3D, (('start', None), ('interior', None), ('end', None)))}

def pipes_from_volmdlr_primitives(primitive, d):
    hy_class, args = VM_EQUIVALENCES[primitive.__class__]
    args2 = []
    for a, index in args:
        if index is None:
            args2.append(getattr(primitive, a))
        else:
            args2.append(getattr(primitive, a)[index])
    return hy_class(*args2, d)

def get_equivalent(pt_coord, known_val, *grid_val):
    """
    Value for pt_coord using linear interpolation from discrete values inside a grid
    y_val can be a list or numpy array
    In case known_val is numpy array dim 2 then grid_val is 2 lists for lines and columns
    """
    dim = len(npy.shape(known_val))
    if dim == 1:
        Lx = grid_val[0]
        x = pt_coord
        y_val = known_val
        n, i = len(Lx), 0
        while i < n and Lx[i] < x:
            i += 1
        if i == n:
            y = y_val[-1]
        else:
            y = y_val[i-1] + (x - Lx[i-1])*(y_val[i] - y_val[i-1])/(Lx[i] - Lx[i-1])
        return y

    elif dim == 2:
        [x, y] = pt_coord
        [Lx, Ly] = grid_val
        nx, i = len(Lx), 0
        while i < nx and Lx[i] < x:
            i += 1
        ny, j = len(Ly), 0
        while j < ny and Ly[i] < y:
            j += 1
        i, j = i-1, j-1
        if i == nx-1 and j == ny-1:
            z = known_val[i, j]
        elif i == nx-1:
            z = known_val[i, j]\
                + (y - Ly[j])*(known_val[i, j+1] - known_val[i, j])/(Ly[j+1] - Ly[j])
        elif j == ny-1:
            z = known_val[i, j]\
                + (x - Lx[i])*(known_val[i+1, j] - known_val[i, j])/(Lx[i+1] - Lx[i])
        else:
            z = known_val[i, j]
            # Interpolation on x
            z += (known_val[i+1, j] - known_val[i, j])*(x - Lx[i])/(Lx[i+1] - Lx[i])
            # Interpolation on y
            z += (known_val[i, j+1] - known_val[i, j])*(y - Ly[j])/(Ly[j+1] - Ly[i])
            # Interpolation on xy
            z += (known_val[i+1, j+1] + known_val[i, j] - known_val[i, j+1] - known_val[i+1, j])\
                *((y - Ly[j])/(Ly[j+1] - Ly[i]))*((x - Lx[i])/(Lx[i+1] - Lx[i]))
        return z
    else:
        print("Function GetEquivalent doesn't support dimension {}".format(dim))
        return None

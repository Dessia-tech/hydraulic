#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 14:24:59 2018

@author: steven
"""

import math
import numpy as npy
import matplotlib.pyplot as plt
import volmdlr as vm
import volmdlr.primitives3D as p3D
from copy import copy

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

class StraightPipe:
    """
    Abstract Class
    Straight pipes linking 2 points
    """
    def __init__(self, p1, p2, d, heat_exchange=True, name=''):
        self.points = [p1, p2]
        self.active_points = self.points
        self.heat_exchange = heat_exchange
        self.radius = d/2
        self.surf = math.pi*self.radius**2
        self.length = p1.PointDistance(p2)
        self.fQ = 16*self.length/(math.pi*self.radius**4)
        self.n_equations = 2
        self.name = name

    def __str__(self):
        return "{} from {} to {}".format(self.__class__.__name__, self.points[0], self.points[1])

    def SystemMatrix(self, constant):
        system_matrix = npy.array([[-constant, self.fQ, constant, 0],
                                   [0, 1, 0, 1]])
        return system_matrix

    def Repr1D(self, j, points_index):
        line = "Line({}) = [{},{}];\n".format(j, points_index[self.points[0]],
                                              points_index[self.points[1]])
        return (line, [j])

    def Dict(self):
        p1, p2 = self.points
        d = {'p1' : p1.Dict(), 'p2' : p2.Dict(), 'd' : self.radius*2,
             'heat_exhcange' : self.heat_exchange, 'length' : float(self.length),
             'fQ' : float(self.fQ), 'name' : self.name}
        return d

class StraightPipe2D(StraightPipe):
    """
    Straight pipes linking 2 2D points
    """
    def __init__(self, p1, p2, d, heat_exchange=True, name=''):
        StraightPipe.__init__(self, p1, p2, d, heat_exchange, name)

    def Draw(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        vm.LineSegment2D(*self.points).MPLPlot(ax)

class StraightPipe3D(StraightPipe):
    """
    Straight pipes linking 2 3D points
    """
    def __init__(self, p1, p2, d, heat_exchange=True, name=''):
        StraightPipe.__init__(self, p1, p2, d, heat_exchange, name)

    def Draw(self, x3D, y3D, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        vm.LineSegment3D(*self.points).MPLPlot2D(x3D, y3D, ax)

    def CADVolume(self):
        axis = self.points[1] - self.points[0]
        axis.Normalize()
        return p3D.HollowCylinder(0.5*(self.points[0] + self.points[1]), axis,
                                  self.radius, self.radius + 0.001,
                                  self.length,
                                  name=self.name)

    @classmethod
    def DictToObject(cls, dict_):
        p1 = vm.Point3D.DictToObject(dict_['p1'])
        p2 = vm.Point3D.DictToObject(dict_['p2'])
        d = dict_['d']
        heat_exchange = dict_['heat_exchange']
        name = dict_['name']
        pipe = cls(p1, p2, d, heat_exchange, name)
        return pipe

class SingularPipe:
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

    def SystemMatrix(self, constant):
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
        length_d = GetEquivalent(abs(self.turn_angle*180/math.pi), ben_LD, ben_angle)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)
        self.n_equations = 2

    def Repr1D(self, j, points_index):
        # Returns the 1D representation of the pipe for gmsh
        n = max(points_index.values())+1
        circle = "Point({0}) = [{1[0]},{1[1]},0.,1.];\n".format(n+1, self.center)
        points_index["center of {}".format(str(self))] = n+1
        circle += "Point({0}) = [{1[0]},{1[1]},0.,1.];\n".format(n+2, self.third_point)
        points_index["intermediate point of {}".format(str(self))] = n+2
        circle += "Circle({}) = [{},{},{}];\n".format(j, points_index[self.points[0]], n+1, n+2)
        circle += "Circle({}) = [{},{},{}];\n".format(j+1, n+2, n+1, points_index[self.points[1]])
        phys = [j, j+1]
        return (circle, phys)

    def Dict(self):
        p1 = self.start_point
        p = self.interior_point
        p2 = self.end_point
        d = {'p1' : p1.Dict(), 'p' : p.Dict(), 'p2' : p2.Dict(),
             'd' : self.radius*2, 'length' : float(self.length),
             'heat_exhcange' : self.heat_exchange,
             'fQ' : float(self.fQ), 'name' : self.name}
        return d

class Bend2D(Bend):
    """
    Bend between two points, defined by third, interior point
    """
    def __init__(self, start_point, interior_point, end_point,
                 diameter,
                 heat_exchange=True,
                 name=''):
        arc = vm.Arc2D(start_point, interior_point, end_point)
        Bend.__init__(self, start_point, interior_point, end_point,
                      arc, diameter, heat_exchange, name)

    def Draw(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        self.arc.MPLPlot(ax)

class Bend3D(Bend):
    """
    Bend between two points, defined by third, interior point
    """
    def __init__(self, start_point, interior_point, end_point,
                 diameter, heat_exchange=True, name=''):
        arc = vm.Arc3D(start_point, interior_point, end_point)
        Bend.__init__(self, start_point, interior_point, end_point,
                      arc, diameter, heat_exchange, name)

    def Draw(self, x3D, y3D, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        self.arc.MPLPlot2D(x3D, y3D, ax)

    def CADVolume(self):
        normal_section = (self.arc.start - self.arc.center).Cross(self.arc.normal)
        section = vm.Contour3D([vm.Circle3D(self.arc.start, self.radius+0.001, normal_section)])
        return p3D.Sweep(section, vm.Wire3D([self.arc]))

    @classmethod
    def DictToObject(cls, dict_):
        p1 = vm.Point3D.DictToObject(dict_['p1'])
        p = vm.Point3D.DictToObject(dict_['p'])
        p2 = vm.Point3D.DictToObject(dict_['p2'])
        d = dict_['d']
        heat_exchange = dict_['heat_exchange']
        name = dict_['name']
        pipe = cls(p1, p, p2, d, heat_exchange, name)
        return pipe

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
        length_d = GetEquivalent(angle*180/math.pi, mbd_LD, mbd_angle)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)

    def Repr1D(self, j, points_index):
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
        length_d = GetEquivalent(d1/d2, enl_LD, enl_rap)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)

    def Repr1D(self, j, points_index):
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
        length_d = GetEquivalent(d1/d2, ctr_LD, ctr_rap)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)

    def Repr1D(self, j, points_index):
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
#    def Repr1D(self, j, points_index):
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
        
    def SystemMatrix(self, constant):
        system_matrix = npy.array([[-constant, self.fQ, constant, 0],
                                   [0, 1, 0, 1]])
        return system_matrix
        
    def CADVolume(self):
        axis = self.points[1] - self.points[0]
        l = axis.Norm()
        axis.Normalize()
        
        return p3D.HollowCylinder(0.5*(self.points[0] + self.points[1]), axis,
                                  0.001, 0.002,
                                  l,
                                  name=self.name)

    def Repr1D(self, j, points_index):
        # Returns the 1D representation of the pipe for gmsh
        line = "Line({}) = [{},{}];\n".format(j, points_index[self.points[0]],
                                              points_index[self.points[1]])
        return (line, [j])
    
    def Draw(self, x3D, y3D, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        vm.LineSegment3D(*self.points).MPLPlot2D(x3D, y3D, ax)

    def Dict(self):
        p1, p2 = self.points
        d = {'p1' : p1.Dict(), 'p2' : p2.Dict(),
             'heat_exhcange' : self.heat_exchange,
             'fQ' : float(self.fQ), 'name' : self.name}
        return d
    
    @classmethod
    def DictToObject(cls, dict_):
        p1 = vm.Point3D.DictToObject(dict_['p1'])
        p2 = vm.Point3D.DictToObject(dict_['p2'])
        fQ = dict_['fQ']
        heat_exchange = dict_['heat_exchange']
        name = dict_['name']
        pipe = cls(p1, p2, fQ, heat_exchange, name)
        return pipe

class JunctionPipe:
    """
    Add pressure drop values
    junction linking 1 pipe to 2+ pipes
    """
    def __init__(self, central_point, other_points, diameter, name=''):
        self.points = [central_point] + other_points
        self.active_points = other_points
        self.central_point = central_point
        self.heat_exchange = False
        self.radius = diameter/2
        self.surf = math.pi*self.radius**2
        self.lengths = [point.PointDistance(central_point) for point in other_points]
        self.n_equations = 1 + len(self.active_points)
        self.name = name

        # To be replaced by values for junction pipes
        length_d = [length/diameter for length in self.lengths]
        self.fQs = [16*2*ld/(math.pi*self.radius**3) for ld in length_d]


    def __str__(self):
        return "{}-jun-{}".format(self.central_point, len(self.active_points))

    def SystemMatrix(self, constant):
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

    def Draw(self, x3D=vm.x3D, y3D=vm.y3D, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        for point in self.active_points:
            vm.LineSegment2D(point, self.central_point).MPLPlot(ax)

    def Repr1D(self, j, points_index):
        # Returns the 1D representation of the pipe for gmsh
        line = "Line({}) = [{},{}];\n".format(j, points_index[self.active_points[0]],
                                              points_index[self.active_points[1]])
        return(line, [j])

VM_EQUIVALENCES = {vm.LineSegment2D: (StraightPipe2D, (('points', 0), ('points', 1))),
                   vm.Arc2D: (Bend2D, (('start', None), ('interior', None), ('end', None))),
                   vm.LineSegment3D: (StraightPipe3D, (('points', 0), ('points', 1))),
                   vm.Arc3D: (Bend3D, (('start', None), ('interior', None), ('end', None)))}

def PipesFromVolmdlrPrimitives(primitive, d):
    hy_class, args = VM_EQUIVALENCES[primitive.__class__]
    args2 = []
    for a, index in args:
        if index is None:
            args2.append(getattr(primitive, a))
        else:
            args2.append(getattr(primitive, a)[index])
    return hy_class(*args2, d)

def GetEquivalent(pt_coord, known_val, *grid_val):
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

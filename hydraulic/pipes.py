#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 14:24:59 2018

@author: steven
"""

import math
import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D


# Definition of equivalent L/D values
# lists in increasing order
ben_angle = [0,45,90,180,360]
ben_LD = [0,16,30,60,120]
# bend radius ?
mbd_angle = [0,45,90]
mbd_LD = [0,15,58]
enl_rap = [1/4,1/2,3/4,1]
enl_LD = [30,20,6.5,0]
ctr_rap = [1,4/3,2,4]
ctr_LD = [0,6.5,11,15]

class StraightPipe:
    """
    Abstract Class
    Straight pipes linking 2 points
    """
    def __init__(self, p1, p2, d, name=''):
        self.points = [p1, p2]
        self.radius = d/2
        self.surf = math.pi*self.radius**2
        self.length = p1.PointDistance(p2)
        self.fQ = 16*self.length/(math.pi*self.radius**4)
        self.name = name
        
    def __repr__(self):
        return("{}-str-{}".format(self.points[0],self.points[1]))
        
    def __str__(self):
        return("{}-str-{}".format(self.points[0],self.points[1]))
        
    def Repr1D(self,j,points_index):
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])


class StraightPipe2D:
    """
    Straight pipes linking 2 2D points
    """
    def __init__(self, p1, p2, d, name=''):
        StraightPipe.__init__(self, p1, p2, d, name)

    def Draw(self, ax):
        l = vm.LineSegment2D(*self.points)
        l.MPLPlot(ax)

        
class StraightPipe3D:
    """
    Straight pipes linking 2 3D points
    """
    def __init__(self, p1, p2, d, name=''):
        StraightPipe.__init__(self, p1, p2, d, name)

    def Draw(self, x3D, y3D, ax):
        l = vm.LineSegment3D(*self.points)
        l.MPLPlot2D(x3D, y3D, ax)    
        
    def CADVolume(self):
        axis = self.points[1]-self.points[0]
        axis.Normalize()
        return primitives3D.HollowCylinder(0.5 * (self.points[0]+self.points[1]), axis,
                                           self.radius, self.radius+0.001,
                                           self.length,
                                           name=self.name)

class SingularPipe:
    """
    other type of pipes linking 2 points
    """
    def __init__(self,p1, p2, form):
        self.points = [p1, p2]
        self.type=form
        
    def __repr__(self):
        return("{}-{}-{}".format(self.points[0],self.type,self.points[1]))
        
    def __str__(self):
        return("{}-{}-{}".format(self.points[0],self.type,self.points[1]))
    
class Bend(SingularPipe):
    """
    Bend between two points, defined by third, interior point
    """
    def __init__(self, start_point, interior_point, end_point, arc, diameter):
        SingularPipe.__init__(self, start_point, end_point, 'ben')
        self.start_point = start_point
        self.interior_point = interior_point
        self.end_point = end_point
        self.radius = 0.5 * diameter
        self.arc = arc
        self.turn_radius = self.arc.radius
        self.turn_angle = self.arc.angle
        self.section = math.pi*self.radius**2
        self.length = self.turn_radius*self.turn_angle
        length_d = GetEquivalent(abs(self.turn_angle*180/math.pi) ,ben_LD, ben_angle)
        self.fQ=16*2*length_d/(math.pi*self.radius**3)
                
        

        
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        n=max(points_index.values())+1
        circle="Point({0}) = [{1[0]},{1[1]},0.,1.];\n".format(n+1,self.center)
        points_index["center of {}".format(str(self))]=n+1
        circle+="Point({0}) = [{1[0]},{1[1]},0.,1.];\n".format(n+2,self.third_point)
        points_index["intermediate point of {}".format(str(self))]=n+2
        circle+="Circle({}) = [{},{},{}];\n".format(j,points_index[self.points[0]],n+1,n+2)
        circle+="Circle({}) = [{},{},{}];\n".format(j+1,n+2,n+1,points_index[self.points[1]])
        phys=[j,j+1]
        return(circle,phys)
        
class Bend2D(Bend):
    """
    Bend between two points, defined by third, interior point
    """
    def __init__(self, start_point, interior_point, end_point, diameter):
        arc = vm.Arc2D(start_point, interior_point, end_point)
        Bend.__init__(self, start_point, interior_point, end_point, arc, diameter)
        
    def Draw(self, ax):
        self.arc.MPLPlot(ax)

class Bend3D(Bend):
    """
    Bend between two points, defined by third, interior point
    """
    def __init__(self, start_point, interior_point, end_point, diameter):
        arc = vm.Arc3D(start_point, interior_point, end_point)
        Bend.__init__(self, start_point, interior_point, end_point, arc, diameter)

    def Draw(self, x3D, y3D, ax):
        self.arc.MPLPlot2D(x3D, y3D, ax)    

    def CADVolume(self):
        normal_section = (self.arc.start-self.arc.center).Cross(self.arc.normal)
        section = vm.Contour3D([vm.Circle3D(self.arc.start, self.radius+0.001, normal_section)])
        return primitives3D.Sweep(section, vm.Wire3D([self.arc]))


class MitterBend(SingularPipe):
    """
    mitter bend 
    (for one point, specify 2 points with same coordinates)
    """
    def __init__(self,p1, p2,diameter,angle):
        SingularPipe.__init__(self,p1, p2,'mbd')
        self.radius=diameter/2
        self.surf=math.pi*self.radius**2
        self.length = p1.PointDistance(p2)
        length_d = GetEquivalent(angle*180/math.pi,mbd_LD,mbd_angle)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)
        
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])
        
class Enlargement(SingularPipe):
    """
    Enlargement from diameter 1 at point 1 to diameter 2 at point 2
    for gradual enlargement, specify the angle
    """
    def __init__(self,p1,p2,d1,d2,angle=None):
        SingularPipe.__init__(self, p1, p2, 'enl')
        self.radius=d1/2
        self.surf=math.pi*(d1/2)**2
        self.length = p1.PointDistance(p2)
        length_d=GetEquivalent(d1/d2,enl_LD,enl_rap)
        self.fQ=16*2*length_d/(math.pi*self.radius**3)
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])
        
class Contraction(SingularPipe):
    """
    Contraction from diameter 1 at point 1 to diameter 2 at point 2
    for gradual ontraction, specify the angle
    """
    def __init__(self,p1,p2,d1,d2,angle=None):
        SingularPipe.__init__(self,p1,p2,'ctr')
        self.radius=d1/2
        self.surf=math.pi*(d1/2)**2
        self.length= p1.PointDistance(p2)
        length_d = GetEquivalent(d1/d2,ctr_LD,ctr_rap)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)
    
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])
        
class UserDefined(SingularPipe):
    """
    Singular Pipe defined by user with a known parameter and different 
    L/D equivalent ratios for different values of the parameter.
    """
    def __init__(self, p1, p2, diameter, param_values, LD_values, *params):
        SingularPipe.__init__(self, p1, p2, 'usr')
        self.radius=diameter/2
        self.surf=math.pi*self.radius**2
        self.length = p1.PointDistance(p2)
        length_d=GetEquivalent(param_values,LD_values,*params)
        self.fQ=16*2*length_d/(math.pi*self.radius**3)
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])
        
class JunctionPipe:
    """
    add pressure drop values
    junction linking 1 pipe to 2+ pipes
    """
    def __init__(self, p1, p2, diameter, value_LD=0):
        self.points = [p1, p2]
        self.radius = diameter/2
        self.surf = math.pi*self.radius**2
        self.length = p1.PointDistance(p2)
        length_d=self.length/diameter # to be replaced by values for junction pipes
        self.fQ=16*2*length_d/(math.pi*self.radius**3)
        
    def __repr__(self):
        return("{}-jun-{}".format(self.points[0],self.points[1]))
        
    def __str__(self):
        return("{}-jun-{}".format(self.points[0],self.points[1]))
        
    def Draw(self, ax):
        vm.LineSegment2D(*self.points).MPLPlot(ax)
        
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])
 
vm_equivalences = {vm.LineSegment3D: (StraightPipe3D, (('points', 0), ('points', 1))),
                   vm.Arc3D: (Bend3D, (('start', None), ('interior', None), ('end', None)))}

def PipesFromVolmdlrPrimitives(primitive, d):
    hy_class, args = vm_equivalences[primitive.__class__]
    args2 = []
    for a, index in args:
        if index is None:
            args2.append(getattr(primitive, a))
        else:
            args2.append(getattr(primitive, a)[index])
    return hy_class(*args2, d)
    
def GetEquivalent(pt_coord, known_val, *grid_val):
    """
    value for pt_coord using linear interpolation from discrete values inside a grid
    y_val can be a list or numpy array
    in case known_val is numpy array dim 2 then grid_val is 2 lists for lines and columns
    """
    dim=len(npy.shape(known_val))
    if dim==1:
        Lx=grid_val[0]
        x=pt_coord
        y_val=known_val
        n,i=len(Lx),0
        while i<n and Lx[i]<x:
            i+=1
        if i==n:
            y=y_val[-1]
        else:
            y=y_val[i-1]+(x-Lx[i-1])*(y_val[i]-y_val[i-1])/(Lx[i]-Lx[i-1])
        return(y)
    elif dim==2:
        [x,y]=pt_coord
        [Lx,Ly]=grid_val
        nx,i=len(Lx),0
        while i<nx and Lx[i]<x:
            i+=1
        ny,j=len(Ly),0
        while j<ny and Ly[i]<y:
            j+=1
        i,j=i-1,j-1
        if i==nx-1 and j==ny-1:
            z=known_val[i,j]
        elif i==nx-1:
            z=known_val[i,j]+(y-Ly[j])*(known_val[i,j+1]-known_val[i,j])/(Ly[j+1]-Ly[j])
        elif j==ny-1:
            z=known_val[i,j]+(x-Lx[i])*(known_val[i+1,j]-known_val[i,j])/(Lx[i+1]-Lx[i])
        else:
            z=known_val[i,j]
            # interpolation on x
            z+=(known_val[i+1,j]-known_val[i,j])*(x-Lx[i])/(Lx[i+1]-Lx[i])
            # interpolation on y
            z+=(known_val[i,j+1]-known_val[i,j])*(y-Ly[j])/(Ly[j+1]-Ly[i])
            # interpolation on xy
            z+=(known_val[i+1,j+1]+known_val[i,j]-known_val[i,j+1]-known_val[i+1,j])*((y-Ly[j])/(Ly[j+1]-Ly[i]))*((x-Lx[i])/(Lx[i+1]-Lx[i]))
        return(z)
    else:
        print("Function GetEquivalent doesn't support dimension {}".format(dim))
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 17:29:11 2018

@author: jezequel
"""
import hydraulic as hy
from hydraulic.fluids import water
import hydraulic.thermal as thermal
import volmdlr as vm
import volmdlr.primitives2d as vmp2d
import volmdlr.primitives3d as vmp3d
import math

diameter = 0.005
dp = 1e5
cp_thickness = 0.001
thermal_conductivity = 237

p0 = vm.Point2D(0, 0)
p1 = vm.Point2D(1, 0)
p2 = vm.Point2D(2, 0)
p3 = vm.Point2D(2, 1)
p4 = vm.Point2D(2, 2)
p5 = vm.Point2D(1, 2)
p6 = vm.Point2D(0, 2)

points = [p0, p1, p2, p3, p4, p5, p6]

rls = vmp2d.OpenedRoundedLineSegments2D(points, {2 : 0.5, 4 : 0.05}, adapt_radius=True)
pipes = [hy.pipes.pipes_from_volmdlr_primitives(p, diameter) for p in rls.basis_primitives]
points = [p.start for p in rls.basis_primitives]+[rls.basis_primitives[-1].end]

bc = [hy.PressureCondition(points[0], dp), hy.PressureCondition(points[-1], 0)]

circuit = hy.Circuit2D(points, pipes, bc, water)
circuit.to_gmsh()

result = circuit.solve_fluidics()
result.circuit.plot()
result.plot()
result.plot(position=True)


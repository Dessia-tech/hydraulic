#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test for importing geometry from volmdlr
"""

import hydraulic as hy
from hydraulic.fluids import water

import volmdlr as vm
import volmdlr.primitives3d as primitives3d

T_liq = 20
dP = 1e5
dQ = 0.005

diameter = 0.010

p1 = vm.Point3D(0, 0, 0)
p2 = vm.Point3D(0.2, 0, 0)
p3 = vm.Point3D(0.2, 0.1, 0)
p4 = vm.Point3D(0.2, 0.1, 0.08)
p5 = vm.Point3D(0.4, 0, 0)

points = [p1, p2, p3, p4, p5]

rl = primitives3d.OpenRoundedLineSegments3D(points, {1: 0.03, 2: 0.01, 3 : 0.02}, adapt_radius=True)

pipes = [hy.pipes.pipes_from_volmdlr_primitives(p, diameter) for p in rl.primitives]
points = [p.points[0] for p in rl.primitives]+[rl.primitives[-1].points[-1]]

# Pipes definitions
bc = [hy.PressureCondition(points[0], dP), hy.PressureCondition(points[-1], 0)]

circuit = hy.Circuit3D(points, pipes, bc, water)

result = circuit.solve_fluidics() 


circuit.plot()
result.plot()

circuit.babylonjs()
# r = circuit.FreeCADExport()
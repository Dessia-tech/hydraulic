#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test for importing geometry from volmdlr
"""

import hydraulic as hy
from hydraulic.fluids import water

import volmdlr as vm
import volmdlr.primitives3D as primitives3D

T_liq = 20
dP = 1e5
dQ = 0.005

diameter = 0.010

p1 = vm.Point3D((0, 0, 0))
p2 = vm.Point3D((0.2, 0, 0))
p3 = vm.Point3D((0.2, 0.1, 0))
p4 = vm.Point3D((0.2, 0.1, 0.08))
p5 = vm.Point3D((0.4, 0, 0))

points = [p1, p2, p3, p4, p5]

rl = primitives3D.RoundedLineSegments3D(points, {1: 0.03, 2: 0.01, 3 : 0.02}, adapt_radius=True)

pipes = [hy.PipesFromVolmdlrPrimitives(p, diameter) for p in rl.basis_primitives]
points = [p.points[0] for p in rl.basis_primitives]+[rl.basis_primitives[-1].points[-1]]

# Pipes definitions

circuit = hy.Circuit3D(points, pipes, water)

pressure, flows = circuit.Solve(imposed = 'pressure', imposed_values = {points[0]: dP, points[-1]: 0}) 

circuit.Draw(vm.x3D, vm.y3D)
circuit.Draw(vm.x3D, vm.z3D)
circuit.DisplaySolution()
#circuit.DisplaySolution(position = 1)

r = circuit.FreeCADExport()
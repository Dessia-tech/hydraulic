#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 17:29:11 2018

@author: jezequel
"""
import hydraulic as hy
from hydraulic.fluids import water
import thermal
import volmdlr as vm
import volmdlr.primitives2D as vmp2D
import volmdlr.primitives3D as vmp3D
import math

diameter = 0.005
dp = 1e5
cp_thickness = 0.001
thermal_conductivity = 237

p0 = vm.Point2D((0, 0, 0))
p1 = vm.Point2D((1, 0, 0))
p2 = vm.Point2D((2, 0, 0))
p3 = vm.Point2D((2, 1, 0))
p4 = vm.Point2D((2, 2, 0))
p5 = vm.Point2D((1, 2, 0))
p6 = vm.Point2D((0, 2, 0))
points = [p0, p1, p2, p3, p4, p5, p6]

rls = vmp2D.RoundedLineSegments2D(points, {2 : 0.5, 4 : 0.05}, adapt_radius=True)
pipes = [hy.pipes.PipesFromVolmdlrPrimitives(p, diameter) for p in rls.basis_primitives]
points = [p.points[0] for p in rls.basis_primitives]+[rls.basis_primitives[-1].points[-1]]

circuit = hy.Circuit2D(points, pipes, water)

result = circuit.SolveFluidics("pressure", {points[0] : dp, points[-1] : 0})
result.circuit.Draw()
result.DisplaySolution()
result.DisplaySolution(position=True)

cooling_plate_parts = []
for i, pipe in enumerate(circuit.pipes):
    thermal_area = 2*math.pi*pipe.radius*pipe.length
    resistor = thermal.ThermalResistor(cp_thickness, thermal_conductivity, thermal_area)
    cooling_plate_parts.append(thermal.CoolingPlatePartFromPipe(pipe, resistor, circuit.fluid))
cooling_plate = thermal.CoolingPlate(circuit, cooling_plate_parts)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 09:34:58 2020

@author: steven
"""

import hydraulic as hy
import hydraulic.thermal as th
from hydraulic.fluids import water
import volmdlr as vm

diameter = 0.010
T_liq = 20
dP = 1e5
dQ = 0.005
constant = 2/(water.rho*water.nu)
coordin = [[0, 0, 0], [.01, 0, 0], [.02, .01, 0], [.02, .03, 0], [.03, .04, 0],
           [.04, .04, 0], [.05, .03, 0], [.05, .01, 0], [.06, 0, 0]]

dl = 0.1

p1 = vm.Point3D((0, 0, 0))
p2 = vm.Point3D((dl, 0.5*dl, 0))
p3 = vm.Point3D((dl, -0.5*dl, 0))
p4 = vm.Point3D((2*dl, 0., 0))
p5 = vm.Point3D((2*dl, 1.1*dl, 0))
p6 = vm.Point3D((3*dl, 0.5*dl, 0.1*dl))
p7 = vm.Point3D((3*dl, -0.5*dl, 0))
p8 = vm.Point3D((2*dl, -1.1*dl, -0.05*dl))


# Middle points
mp123 = (p1+p2+p3)/3.
mp245 = (p2+p4+p5)/3.
mp348 = (p3+p4+p8)/3.
mp467 = (p4+p6+p7)/3.

points = [p1, p2, p3, p4, p5, p6, p7, p8, mp123, mp245, mp348, mp467]

# Pipes definitions
junction1 = hy.pipes.JunctionPipe(mp123, [p1, p2, p3], diameter)
junction2 = hy.pipes.JunctionPipe(mp245, [p2, p4, p5], diameter)
junction3 = hy.pipes.JunctionPipe(mp348, [p3, p4, p8], diameter)
junction4 = hy.pipes.JunctionPipe(mp467, [p4, p6, p7], diameter)
straight1 = hy.pipes.StraightPipe3D(p5, p6, diameter)
straight2 = hy.pipes.StraightPipe3D(p7, p8, diameter)
pipes = [junction1, junction2, junction3, junction4, straight1, straight2]

# pipes = [hy.pipes.StraightPipe3D(points[link[0]], points[link[1]], diameter, False) for link in l_strt_link]


boundary_conditions = [hy.PressureCondition(p1, dP), hy.PressureCondition(p4, 0)]
circuit = hy.Circuit3D(points, pipes, boundary_conditions, water)
circuit.plot()
circuit.plot_graph()

volume_model = circuit.volmdlr_volume_model()
volume_model.babylonjs()

# Fluidics calculations
fluidics_result = circuit.SolveFluidics()
fluidics_result.DisplaySolution()
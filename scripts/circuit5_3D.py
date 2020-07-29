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
l_coord = [[0, 0, 0], [.01, 0, 0], [.02, .01, 0], [.02, .03, 0], [.03, .04, 0],
           [.04, .04, 0], [.05, .03, 0], [.05, .01, 0], [.06, 0, 0]]
points = [vm.Point3D(coord) for coord in l_coord]
l_strt_link = [[4, 7], [6, 8], [8, 5]]
l_junc_link = {0: [1, 2, 3], 1: [0, 2, 4], 2: [0, 1, 6], 3: [0, 7, 5]}
# Pipes definitions
pipes = [hy.pipes.StraightPipe3D(points[link[0]], points[link[1]], diameter, False) for link in l_strt_link]
for central_point_index, other_points_indices in l_junc_link.items():
    other_points = [points[i] for i in other_points_indices]
    pipes.append(hy.pipes.JunctionPipe(points[central_point_index],
                                       other_points,
                                       diameter))
boundary_conditions = [hy.PressureCondition(points[0], dP), hy.PressureCondition(points[3], 0)]
circuit = hy.Circuit3D(points, pipes, boundary_conditions, water)
# Fluidics calculations
fluidics_result = circuit.SolveFluidics()
fluidics_result.DisplaySolution()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 12:24:00 2018

@author: steven
"""

import hydraulic as hy
from hydraulic.fluids import water
import volmdlr as vm

diameter = 0.008

T_liq = 20
dP = 1e5
dQ = 0.005

l_coord = [[0,0],[.01,0],[.02,.01],[.02,.03],[.03,.04],[.04,.04],[.05,.03],
         [.05,.01],[.06,0],[.06,.01],[.06,.03],[.07,.04],[.08,.04],
         [.09,.04],[.10,.04],[.07,0],[.07,.01],[.07,.02],[.08,.03],[.09,.03],
         [.08,0],[.10,0],[.10,.01],[.09,.02],[.09,0],[.11,.03],[.11,0],
         [.10,-0.01],[.09,-.01],[.08,-.01],[0,-.01]]
points = [vm.Point2D(coord) for coord in l_coord]

l_strt_link = [[0,1],[2,3],[4,5],[6,7],[9,10],[11,12],[13,14],[16,17],[24,21],
               [25,26],[27,28],[29,30],[28,29],[20,24]]
l_junc_link = [[8,15],[8,9],[12,13],[19,13],[15,16],[15,20],[18,19],[23,19]]

# Pipes definitions
pipes = [hy.pipes.StraightPipe2D(points[link[0]], points[link[1]], diameter) for link in l_strt_link]
pipes.append(hy.pipes.Bend2D(points[1], vm.Point2D((0.017,0.003)), points[2], diameter))
pipes.append(hy.pipes.Bend2D(points[3], vm.Point2D((0.023,0.037)), points[4], diameter))
pipes.append(hy.pipes.Bend2D(points[5], vm.Point2D((0.047,0.037)), points[6], diameter))
pipes.append(hy.pipes.Bend2D(points[7], vm.Point2D((0.053,0.003)), points[8], diameter))
pipes.append(hy.pipes.Bend2D(points[10], vm.Point2D((0.063,0.037)), points[11], diameter))
pipes.append(hy.pipes.Bend2D(points[17], vm.Point2D((0.073,0.027)), points[18], diameter))
pipes.append(hy.pipes.Bend2D(points[21], vm.Point2D((0.105,0.005)), points[22], diameter))
pipes.append(hy.pipes.Bend2D(points[22], vm.Point2D((0.093,0.013)), points[23], diameter))
pipes.append(hy.pipes.Bend2D(points[14], vm.Point2D((0.107,0.037)), points[25], diameter))
pipes.append(hy.pipes.Bend2D(points[26], vm.Point2D((0.107,-0.007)), points[27], diameter))
pipes.extend([hy.pipes.JunctionPipe(points[link[0]], points[link[1]], diameter) for link in l_junc_link])

delta_P = {points[30]:0, points[0]: dP}
Qs = {points[30]: "out", points[0]: dQ}

circuit = hy.Circuit2D(points, pipes, water)

bords = [[0,-0.05],[0.12,-0.05],[0.12,0.05],[0,0.05]]
liq_T = {points[0]:T_liq}

result = circuit.SolveFluidics(imposed = 'pressure', imposed_values = {points[0]: dP, points[-1]: 0}) 

result.circuit.Draw()
result.DisplaySolution()
result.DisplaySolution(position=True)
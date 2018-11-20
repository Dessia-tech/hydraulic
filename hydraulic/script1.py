#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 17:24:47 2018

@author: jezequel
"""

import thermal as th
from hydraulic.fluids import water
import math

# =============================================================================
# Resistors in series
# =============================================================================
#node1 = th.Node()
#node2 = th.Node()
#node3 = th.Node()
#nodes = [node1, node2, node3]
#
#bounds1 = th.HeatFlowInBound([node1], -10)
#block1 = th.Resistor([node1, node2], 2, 1, 1)
#block2 = th.Resistor([node2, node3], 10, 1, 1)
#bounds4 = th.TemperatureBound([node3], 313)
#blocks = [bounds1, block1, block2, bounds4]
#
#thc = th.Circuit(nodes, blocks)
#vd = thc.NumberElements()
#junction_input_flows = {}
#system_matrix = thc.SystemMatrix(junction_input_flows)
#system_matrix, vector_b, vect_solution = thc.Solve(junction_input_flows)

# =============================================================================
# Resistors : star disposition
# =============================================================================
#node1 = th.Node()
#node2 = th.Node()
#node3 = th.Node()
#node4 = th.Node()
#nodes = [node1, node2, node3, node4]
#
#bounds1 = th.TemperatureBound([node1], 273)
#block1 = th.Resistor([node1, node2], 2, 1, 1)
#block2 = th.Resistor([node2, node3], 10, 1, 1)
#block3 = th.Resistor([node2, node4], 5, 1, 1)
#bounds4 = th.TemperatureBound([node3], 313)
#bounds5 = th.TemperatureBound([node4], 313)
#blocks = [bounds1, block1, block2, block3, bounds4, bounds5]
#
#thc = th.Circuit(nodes, blocks)
#vd = thc.NumberElements()
#junction_input_flows = {}
#system_matrix = thc.SystemMatrix(junction_input_flows)
#system_matrix, vector_b, vect_solution = thc.Solve(junction_input_flows)

# =============================================================================
# Resistors in parallel
# =============================================================================
#node1 = th.Node()
#node2 = th.Node()
#nodes = [node1, node2]
#
#bounds1 = th.HeatFlowInBound([node1], 680)
#block1 = th.Resistor([node1, node2], 2, 1, 1)
#block2 = th.Resistor([node1, node2], 10, 1, 1)
#block3 = th.Resistor([node1, node2], 5, 1, 1)
#bounds4 = th.TemperatureBound([node2], 313)
#blocks = [bounds1, block1, block2, block3, bounds4]
#
#thc = th.Circuit(nodes, blocks)
#vd = thc.NumberElements()
#junction_input_flows = {}
#system_matrix = thc.SystemMatrix(junction_input_flows)
#system_matrix, vector_b, vect_solution = thc.Solve(junction_input_flows)

# =============================================================================
# !!! Resistors in square
# =============================================================================
#node0 = th.Node()
#node1 = th.Node()
#node2 = th.Node()
#node3 = th.Node()
#node4 = th.Node()
#nodes = [node0, node1, node2, node3, node4]
#
#bounds1 = th.HeatFlowInBound([node0], 0)
#block1 = th.Resistor([node0, node1], 2, 1, 1)
#block2 = th.Resistor([node0, node2], 10, 1, 1)
#block3 = th.Resistor([node0, node3], 5, 1, 1)
#block4 = th.Resistor([node0, node4], 1, 1, 1)
#bounds5 = th.HeatFlowOutBound([node1])
#bounds6 = th.HeatFlowOutBound([node2])
#bounds7 = th.HeatFlowOutBound([node3])
#bounds8 = th.HeatFlowOutBound([node4])
#blocks = [bounds1, block1, block2, block3, block4, bounds5, bounds6, bounds7, bounds8]
#
#thc = th.Circuit(nodes, blocks)
#vd = thc.NumberElements()
#junction_input_flows = {}
#system_matrix = thc.SystemMatrix(junction_input_flows)
#system_matrix, vector_b, vect_solution = thc.Solve(junction_input_flows)

# =============================================================================
# Cooling pipe
# =============================================================================
#node1 = th.Node()
#node2 = th.Node()
#node3 = th.Node()
#node4 = th.Node()
#nodes = [node1, node2, node3, node4]
#
#bounds1 = th.TemperatureBound([node1], 273)
#block1 = th.ThermalPipe([node1, node2, node3])
#bounds2 = th.HeatFlowOutBound([node2])
#block2 = th.Resistor([node3, node4], 2, 1, 1)
#bounds3 = th.HeatFlowInBound([node4], -100)
#blocks = [bounds1, block1, bounds2, block2, bounds3]
#
#thc = th.Circuit(nodes, blocks)
#vd = thc.NumberElements()
#input_flows = {block1 : [0.000001]}
#system_matrix = thc.SystemMatrix(input_flows, water)
#system_matrix, vector_b, vect_solution = thc.Solve(input_flows, water)

# =============================================================================
# Double Cooling Pipe
# =============================================================================
node01 = th.Node()
node02 = th.Node()
node03 = th.Node()
node04 = th.Node()
node11 = th.Node()
node12 = th.Node()
node13 = th.Node()
node14 = th.Node()
nodes = [node01, node02, node03, node04, node11, node12, node13, node14]

bounds01 = th.TemperatureBound([node01], 273)
block01 = th.ThermalPipe([node01, node02, node03])
block02 = th.Resistor([node03, node04], 2, 1, 1)
bounds02 = th.HeatFlowInBound([node04], -100)

link0 = th.NodeEquivalence([node02, node11])

block11 = th.ThermalPipe([node11, node12, node13])
bounds12 = th.HeatFlowOutBound([node12])
block12 = th.Resistor([node13, node14], 5, 1, 1)
bounds13 = th.HeatFlowInBound([node14], -200)
blocks = [bounds01, block01, block02, bounds02, link0 , block11, bounds12, block12, bounds13]

thc = th.Circuit(nodes, blocks)
input_flows = {block01 : [0.000001], block11 : [0.000001]}
system_matrix = thc.SystemMatrix(input_flows, water)
system_matrix, vector_b, vect_solution = thc.Solve(input_flows, water)

# =============================================================================
# Double Cooling Pipe variant (without NodeEquivalence)
# =============================================================================
node01 = th.Node()
node02 = th.Node()
node03 = th.Node()
node04 = th.Node()
node12 = th.Node()
node13 = th.Node()
node14 = th.Node()
nodes = [node01, node02, node03, node04, node12, node13, node14]

bounds01 = th.TemperatureBound([node01], 273)
block01 = th.ThermalPipe([node01, node02, node03])
block02 = th.Resistor([node03, node04], 2, 1, 1)
bounds02 = th.HeatFlowInBound([node04], -100)

block11 = th.ThermalPipe([node02, node12, node13])
bounds12 = th.HeatFlowOutBound([node12])
block12 = th.Resistor([node13, node14], 5, 1, 1)
bounds13 = th.HeatFlowInBound([node14], -200)
blocks = [bounds01, block01, block02, bounds02, block11, bounds12, block12, bounds13]

thc = th.Circuit(nodes, blocks)
input_flows = {block01 : [0.000001], block11 : [0.000001]}
system_matrix1 = thc.SystemMatrix(input_flows, water)
system_matrix1, vector_b1, vect_solution1 = thc.Solve(input_flows, water)

# =============================================================================
# Simple Junction
# =============================================================================
#node1 = th.Node()
#node2 = th.Node()
#node3 = th.Node()
#node4 = th.Node()
#node5 = th.Node()
#nodes = [node1, node2, node3, node4, node5]
#
#bounds1 = th.TemperatureBound([node1], 313)
#bounds2 = th.TemperatureBound([node2], 353)
#block = th.Junction([node1, node2], [node3, node4, node5])
#bounds3 = th.HeatFlowOutBound([node3])
#bounds4 = th.HeatFlowOutBound([node4])
#bounds5 = th.HeatFlowOutBound([node5])
#blocks = [bounds1, bounds2, block, bounds3, bounds4, bounds5]
#
#thc = th.Circuit(nodes, blocks)
#vd = thc.NumberElements()
#input_flows = {block : [0.000001, 0.000005]}
#system_matrix = thc.SystemMatrix(input_flows, water)
#system_matrix, vector_b, vect_solution = thc.Solve(input_flows, water)

# =============================================================================
# Multiple Junctions and Pipes
# =============================================================================
#node1 = th.Node()
#node2 = th.Node()
#node3 = th.Node()
#node4 = th.Node()
#node5 = th.Node()
#nodes = [node1, node2, node3, node4, node5]
#
#bounds1 = th.TemperatureBound([node1], 313)
#bounds2 = th.TemperatureBound([node2], 353)
#block = th.Junction([node1, node2], [node3, node4, node5])
#bounds3 = th.HeatFlowOutBound([node3])
#bounds4 = th.HeatFlowOutBound([node4])
#bounds5 = th.HeatFlowOutBound([node5])
#blocks = [bounds1, bounds2, block, bounds3, bounds4, bounds5]
#
#thc = th.Circuit(nodes, blocks)
#vd = thc.NumberElements()
#input_flows = {block : [0.000001, 0.000005]}
#system_matrix = thc.SystemMatrix(input_flows, water)
#system_matrix, vector_b, vect_solution = thc.Solve(input_flows, water)

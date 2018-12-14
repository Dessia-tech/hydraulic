#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 17:24:47 2018

@author: jezequel
"""

import thermal as th
from hydraulic.fluids import water
import math
import volmdlr as vm

# =============================================================================
# Environment
# =============================================================================
cp_thickness = 0.001
thermal_conductivity = 237

# =============================================================================
# Resistors in series
# =============================================================================
#node1 = th.Node()
#node2 = th.Node()
#node3 = th.Node()
#nodes = [node1, node2, node3]
#
#bound1 = th.HeatFlowInBound([node1], -10)
#block1 = th.Resistor([node1, node2], 2, 1, 1)
#block2 = th.Resistor([node2, node3], 10, 1, 1)
#bound4 = th.TemperatureBound([node3], 313)
#blocks = [bound1, block1, block2, bound4]
#
#thc = th.Circuit(nodes, blocks)
#junction_input_flows = {}
#system_matrix = thc.SystemMatrix(junction_input_flows)
#result = thc.Solve(junction_input_flows)

# =============================================================================
# Resistors : star disposition
# =============================================================================
#node1 = th.Node()
#node2 = th.Node()
#node3 = th.Node()
#node4 = th.Node()
#nodes = [node1, node2, node3, node4]
#
#bound1 = th.TemperatureBound([node1], 273)
#block1 = th.Resistor([node1, node2], 2, 1, 1)
#block2 = th.Resistor([node2, node3], 10, 1, 1)
#block3 = th.Resistor([node2, node4], 5, 1, 1)
#bound4 = th.TemperatureBound([node3], 313)
#bound5 = th.TemperatureBound([node4], 313)
#blocks = [bound1, block1, block2, block3, bound4, bound5]
#
#thc = th.Circuit(nodes, blocks)
#junction_input_flows = {}
#system_matrix = thc.SystemMatrix(junction_input_flows)
#result = thc.Solve(junction_input_flows)

# =============================================================================
# Resistors in parallel
# =============================================================================
#node1 = th.Node()
#node2 = th.Node()
#nodes = [node1, node2]
#
#bound1 = th.HeatFlowInBound([node1], 680)
#block1 = th.Resistor([node1, node2], 2, 1, 1)
#block2 = th.Resistor([node1, node2], 10, 1, 1)
#block3 = th.Resistor([node1, node2], 5, 1, 1)
#bound4 = th.TemperatureBound([node2], 313)
#blocks = [bound1, block1, block2, block3, bound4]
#
#thc = th.Circuit(nodes, blocks)
#junction_input_flows = {}
#system_matrix = thc.SystemMatrix(junction_input_flows)
#result = thc.Solve(junction_input_flows)

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
#bound1 = th.HeatFlowInBound([node0], 0)
#block1 = th.Resistor([node0, node1], 2, 1, 1)
#block2 = th.Resistor([node0, node2], 10, 1, 1)
#block3 = th.Resistor([node0, node3], 5, 1, 1)
#block4 = th.Resistor([node0, node4], 1, 1, 1)
#bound5 = th.HeatFlowOutBound([node1])
#bound6 = th.HeatFlowOutBound([node2])
#bound7 = th.HeatFlowOutBound([node3])
#bound8 = th.HeatFlowOutBound([node4])
#blocks = [bound1, block1, block2, block3, block4, bound5, bound6, bound7, bound8]
#
#thc = th.Circuit(nodes, blocks)
#junction_input_flows = {}
#system_matrix = thc.SystemMatrix(junction_input_flows)
#result = thc.Solve(junction_input_flows)

# =============================================================================
# Test on Node equivalence
# =============================================================================
#node1 = th.Node()
#node2 = th.Node()
#node3 = th.Node()
#node4 = th.Node()
#node5 = th.Node()
#node6 = th.Node()
#nodes = [node1, node2, node3, node5, node6]
#
#bound1 = th.TemperatureBound([node1], 273)
#block1 = th.ThermalPipe([node1, node2, node3], 0.0001)
#bound2 = th.HeatFlowOutBound([node3])
##ne1 = th.NodeEquivalence([node2, node4])
#block2 = th.ThermalPipe([node2, node5, node6], 0.0001)
#bound2 = th.HeatFlowOutBound([node6])
#bound2 = th.HeatFlowOutBound([node5])
#blocks = [bound1, block1, block2, bound2]
#
#thc = th.Circuit(nodes, blocks)
#system_matrix = thc.SystemMatrix()
# =============================================================================
# Cooling pipe
# =============================================================================
#node1 = th.Node()
#node2 = th.Node()
#node3 = th.Node()
#node4 = th.Node()
#nodes = [node1, node2, node3, node4]
#
#bound1 = th.TemperatureBound([node1], 273)
#block1 = th.ThermalPipe([node1, node2, node3])
#bound2 = th.HeatFlowOutBound([node2])
#block2 = th.Resistor([node3, node4], 2, 1, 1)
#bound3 = th.HeatFlowInBound([node4], -100)
#blocks = [bound1, block1, bound2, block2, bound3]
#
#thc = th.Circuit(nodes, blocks)
#input_flows = {block1 : [0.000001]}
#system_matrix = thc.SystemMatrix(input_flows, water)
#result = thc.Solve(input_flows, water)

# =============================================================================
# Double Cooling Pipe
# =============================================================================
#node01 = th.Node()
#node02 = th.Node()
#node03 = th.Node()
#node04 = th.Node()
#node11 = th.Node()
#node12 = th.Node()
#node13 = th.Node()
#node14 = th.Node()
#nodes = [node01, node02, node03, node04, node11, node12, node13, node14]
#
#bound01 = th.TemperatureBound([node01], 273)
#block01 = th.ThermalPipe([node01, node02, node03], 0.000001, water)
#block02 = th.Resistor([node03, node04], 2, 1, 1)
#bound02 = th.HeatFlowInBound([node04], -100)
#
#link0 = th.NodeEquivalence([node02, node11])
#
#block11 = th.ThermalPipe([node11, node12, node13], 0.000001, water)
#bound12 = th.HeatFlowOutBound([node12])
#block12 = th.Resistor([node13, node14], 5, 1, 1)
#bound13 = th.HeatFlowInBound([node14], -200)
#blocks = [bound01, block01, block02, bound02, link0 , block11, bound12, block12, bound13]
#
#thc = th.Circuit(nodes, blocks)
#system_matrix = thc.SystemMatrix()
#result = thc.Solve()

# =============================================================================
# Double Cooling Pipe variant (without NodeEquivalence)
# =============================================================================
#node01 = th.Node('h01')
#node02 = th.Node('h02')
#node03 = th.Node('t01')
#node04 = th.Node('t02')
#node12 = th.Node('h12')
#node13 = th.Node('t11')
#node14 = th.Node('t12')
#nodes = [node01, node02, node03, node04, node12, node13, node14]
#
#bound01 = th.TemperatureBound([node01], 273, 'bc0')
#block01 = th.ThermalPipe([node01, node02, node03], 0.000001, water, 'p0')
#block02 = th.Resistor([node03, node04], 2, 1, 'r0')
#bound02 = th.HeatFlowInBound([node04], -100, 'cell0')
#
#block11 = th.ThermalPipe([node02, node12, node13], 0.000001, water, 'p1')
#bound12 = th.HeatFlowOutBound([node12], 'bc1')
#block12 = th.Resistor([node13, node14], 5, 1, 'r1')
#bound13 = th.HeatFlowInBound([node14], -200, 'cell1')
#blocks = [bound01, block01, block02, bound02, block11, bound12, block12, bound13]
#
#thc = th.Circuit(nodes, blocks)
#system_matrix = thc.SystemMatrix()
#result = thc.Solve()
#result.Display()

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
#bound1 = th.TemperatureBound([node1], 313)
#bound2 = th.TemperatureBound([node2], 353)
#block = th.Junction([node1, node2], [node3, node4, node5], nodes, [0.000001, 0.000005], water)
#bound3 = th.HeatFlowOutBound([node3])
#bound4 = th.HeatFlowOutBound([node4])
#bound5 = th.HeatFlowOutBound([node5])
#blocks = [bound1, bound2, block, bound3, bound4, bound5]
#
#thc = th.Circuit(nodes, blocks)
#system_matrix = thc.SystemMatrix()
#result = thc.Solve()

# =============================================================================
# Multiple Junctions and Pipes
# =============================================================================
# Hydraulic
hnodes = [th.Node('h'+str(i)) for i in range(22)]
tnodes = [th.Node('t'+str(i)) for i in range(15)]

b1 = th.TemperatureBound([hnodes[0]], 313, 'b1')
p1 = th.ThermalPipe([hnodes[0], hnodes[1], tnodes[0]], 0.0000030, water, 'p1')
j1 = th.Junction([hnodes[1]], [hnodes[2], hnodes[3]], [0.0000030], water, 'j1')
p2 = th.ThermalPipe([hnodes[2], hnodes[4], tnodes[1]], 0.0000010, water, 'p2')
p3 = th.ThermalPipe([hnodes[4], hnodes[6], tnodes[2]], 0.0000020, water, 'p3')
p4 = th.ThermalPipe([hnodes[6], hnodes[8], tnodes[3]], 0.0000010, water, 'p4')

p5 = th.ThermalPipe([hnodes[3], hnodes[5], tnodes[4]], 0.0000020, water, 'p5')
p6 = th.ThermalPipe([hnodes[5], hnodes[7], tnodes[5]], 0.0000010, water, 'p6')
p7 = th.ThermalPipe([hnodes[7], hnodes[9], tnodes[6]], 0.0000020, water, 'p7')
j2 = th.Junction([hnodes[8], hnodes[9]], [hnodes[10], hnodes[11], hnodes[12]],
                 [0.0000015, 0.0000015], water, 'j2')
p8 = th.ThermalPipe([hnodes[10], hnodes[13], tnodes[7]], 0.0000010, water, 'p8')
p9 = th.ThermalPipe([hnodes[13], hnodes[15], tnodes[8]], 0.0000010, water, 'p9')
p10 = th.ThermalPipe([hnodes[15], hnodes[17], tnodes[9]], 0.0000010, water, 'p10')

p11 = th.ThermalPipe([hnodes[11], hnodes[18], tnodes[10]], 0.0000010, water, 'p11')

p12 = th.ThermalPipe([hnodes[12], hnodes[14], tnodes[11]], 0.0000010, water, 'p12')
p13 = th.ThermalPipe([hnodes[14], hnodes[16], tnodes[12]], 0.0000010, water, 'p13')
p14 = th.ThermalPipe([hnodes[16], hnodes[19], tnodes[13]], 0.0000010, water, 'p14')
j3 = th.Junction([hnodes[17], hnodes[18], hnodes[19]], [hnodes[20]],
                 [0.0000010, 0.0000010, 0.0000010],
                 water, 'j2')
p15 = th.ThermalPipe([hnodes[20], hnodes[21], tnodes[14]], 0.0000030, water, 'p15')
b2 = th.HeatFlowOutBound([hnodes[21]], 'b2')

bounds = [b1, b2]
pipes = [p1, p2, p3, p4, p5, p6, p7,
         p8, p9, p10, p11, p12, p13, p14,
         p15]
junctions = [j1, j2, j3]

hblocks = [*bounds, *pipes, *junctions]

# Thermal
rnodes = []
tblocks = []
for i, pipe in enumerate(pipes):
    tnode = tnodes[i]
    rnode = th.Node('r'+str(i))
    rnodes.append(rnode)
    r = th.Resistor([tnode, rnode], 0.001/0.01, 0.010)
    tblocks.append(r)
    if i > 0 and i <= 3:
        hf = th.HeatFlowInBound([rnode], 10)
    else:
        hf = th.HeatFlowInBound([rnode], -10)
    tblocks.append(hf)

# Circuit
nodes = hnodes + tnodes + rnodes
blocks = hblocks + tblocks
thc = th.Circuit(nodes, blocks)
system_matrix = thc.SystemMatrix()

# Solve and diplay
result = thc.Solve()
result.Display()    

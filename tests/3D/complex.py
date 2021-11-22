#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 14:56:55 2018

@author: jezequel
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

coordinates = [[0, 0, 0], [.01, 0, 0], [.02, .01, 0], [.02, .03, 0], [.03, .04, 0],
               [.04, .04, 0], [.05, .03, 0], [.05, .01, 0], [.06, 0, 0], [.07, 0, 0],
               [.075, .0, .0], [.07, .01, 0], [.07, .03, 0], [.08, .04, 0], [.09, .04, 0],
               [.10, .035, 0], [.10, .04, 0], [.11, .04, 0], [.08, 0, 0], [.08, .01, 0],
               [.08, .02, 0], [.09, .03, 0], [.10, .03, 0], [.09, 0, 0], [.11, 0, 0],
               [.11, .01, 0], [.10, .02, 0], [.10, 0, 0], [.12, .03, 0], [.12, 0, 0],
               [.11, -0.01, 0], [.10, -.01, 0], [.09, -.01, 0], [0, -.01, 0]]
points = [vm.Point3D(*coord) for coord in coordinates]

straight_links = [[0, 1], [2, 3], [4, 5], [6, 7], [11, 12], [13, 14], [19, 20],
               [27, 24], [28, 29], [30, 31], [32, 33], [31, 32], [23, 27]]
junction_links = {9 : [8, 10, 11], 16 : [14, 17, 15], 18 : [10, 19, 23], 22 : [15, 21, 26]}

# Pipes definitions
pipes = [hy.pipes.StraightPipe3D(points[link[0]], points[link[1]], diameter, True) for link in straight_links]
pipes.append(hy.pipes.Bend3D(points[1], vm.Point3D(0.017,0.003, 0), points[2], diameter, False))
pipes.append(hy.pipes.Bend3D(points[3], vm.Point3D(0.023,0.037, 0), points[4], diameter, False))
pipes.append(hy.pipes.Bend3D(points[5], vm.Point3D(0.047,0.037, 0), points[6], diameter, False))
pipes.append(hy.pipes.Bend3D(points[7], vm.Point3D(0.053,0.003, 0), points[8], diameter, False))
pipes.append(hy.pipes.Bend3D(points[12], vm.Point3D(0.073,0.037, 0), points[13], diameter, False))
pipes.append(hy.pipes.Bend3D(points[20], vm.Point3D(0.083,0.027, 0), points[21], diameter, False))
pipes.append(hy.pipes.Bend3D(points[24], vm.Point3D(0.115,0.005, 0), points[25], diameter, False))
pipes.append(hy.pipes.Bend3D(points[25], vm.Point3D(0.103,0.013, 0), points[26], diameter, False))
pipes.append(hy.pipes.Bend3D(points[17], vm.Point3D(0.117,0.037, 0), points[28], diameter, False))
pipes.append(hy.pipes.Bend3D(points[29], vm.Point3D(0.117,-0.007, 0), points[30], diameter, False))

for central_point_index, other_points_indices in junction_links.items():
    other_points = [points[i] for i in other_points_indices]
    pipes.append(hy.pipes.JunctionPipe(points[central_point_index],
                                       other_points,
                                       diameter))
boundary_conditions = [hy.PressureCondition(points[0], dP), hy.PressureCondition(points[-1], 0)]

circuit = hy.Circuit3D(points, pipes, boundary_conditions, water)
circuit.plot()
circuit.plot2d()


# Fluidics calculations
fluidics_result = circuit.solve_fluidics()

# Hydraulic to thermal
thermohydraulic_circuit = fluidics_result.to_thermal([points[0]], [points[-1]])

# Add Resistors
for i, (pipe, block) in enumerate(thermohydraulic_circuit.pipe2block.items()):
    wall_node = block.nodes[-1]
    if pipe.heat_exchange\
    and wall_node in thermohydraulic_circuit.interface_nodes['wall_nodes']:
        condition_node = th.Node('thf'+str(i))
        hf_condition = th.HeatFlowInBound([condition_node], -100, 'bc'+str(i))
        thermohydraulic_circuit.thermal_circuit.add_resistor(wall_node, 10, 1, hf_condition)

# Add boundary conditions
boundary_conditions = []
for input_node in thermohydraulic_circuit.interface_nodes['input']:
    temperature_condition = th.TemperatureBound([input_node], 293, 'temperature_condition')
    boundary_conditions.append(temperature_condition)

for output_node in thermohydraulic_circuit.interface_nodes['output']:
    hf_condition = th.HeatFlowOutBound([output_node], 'hf_condition')
    boundary_conditions.append(hf_condition)
    
thermohydraulic_circuit.thermal_circuit.add_blocks(boundary_conditions)
thermal_results = thermohydraulic_circuit.thermal_circuit.solve()
thermal_results.plot()

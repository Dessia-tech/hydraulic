#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 14:55:10 2019

@author: steven
"""

kr1 = 0.1
kr2 = 0.1
l = 0.25
kc = 10.
S = 0.01
q = -100.

import hydraulic.thermal as thermal

node_fluid_left = thermal.Node('Fluid left')
node_surface_left = thermal.Node('Surface left')
node_surface_right = thermal.Node('Fluid right')
node_fluid_right = thermal.Node('Surface right')
node_inner_medium = thermal.Node('Inner medium')

nodes = [node_fluid_left, node_fluid_right, node_surface_left, node_surface_right, node_inner_medium]

res1 = thermal.Resistor([node_fluid_left, node_surface_left], kr1, 1)
res2 = thermal.Resistor([node_fluid_right, node_surface_right], kr2, 1)
medium = thermal.UnidimensionalMedium([node_surface_left, node_surface_right, node_inner_medium], kc, l, S, 'medium')
imp_temp1 = thermal.TemperatureBound([node_fluid_left], 40)
imp_temp2 = thermal.TemperatureBound([node_fluid_right], 20)
imp_flux = thermal.HeatFlowInBound([node_inner_medium], q)

blocks = [res1, res2, medium, imp_temp1, imp_temp2, imp_flux]
circuit = thermal.Circuit(nodes, blocks)
result = circuit.Solve()
result.Display()
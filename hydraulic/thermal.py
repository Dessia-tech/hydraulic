#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:54:15 2018

@author: jezequel
"""
import numpy as npy
import networkx as nx
from scipy.linalg import solve
from hydraulic.fluids import water

class Node:
    """
    Defines a thermal node
    """
    def __init__(self):
        self.n_equations = 1

class Resistor:
    """
    Defines a thermal resistor
    """
    def __init__(self, nodes, thickness, thermal_conductivity, area):
        self.nodes = nodes
        self.thickness = thickness
        self.thermal_conductivity = thermal_conductivity
        self.area = area
        self.resistance = thickness/(thermal_conductivity*area)
        self.system_matrix = self.SystemMatrix()
        self.n_equations, self.n_variables = self.system_matrix.shape

    def SystemMatrix(self):
        """
        Computes node equations system
        sum(phi) = 0
        """
        matrix = npy.array([[-self.resistance, 1, self.resistance, 0],
                            [0, 1, 0, 1]])
        return matrix

class Junction:
    """
    Defines a junction, with several inputs and outputs
    """
    def __init__(self, input_nodes, output_nodes):
        self.input_nodes = input_nodes
        self.output_nodes = output_nodes
        self.nodes = input_nodes + output_nodes
        input_flows = [1 for i in range(len(input_nodes))]
        self.system_matrix = self.SystemMatrix(input_flows)
        self.n_equations, self.n_variables = self.system_matrix.shape

    def SystemMatrix(self, input_flows):
        """
        Computes block equations system
        Output temperatures are a weighted mean of input volumetric flows
        Heat flows are equals
        """
        matrix_generator = []
        line_generator = []

        # Line 1 : Weigthed mean volumetric flows equation
        for n_node in range(len(self.input_nodes)):
            flow = input_flows[n_node]
            line_generator.extend([flow/sum(input_flows), 0])

        output_generator = [0]*(len(self.output_nodes)*2)
        output_generator[0] = -1
        line_generator.extend(output_generator)
        matrix_generator.append(line_generator)

        # Output temperatures equalities
        for n_node in range(len(self.output_nodes)):
            line_generator = [0]*(len(self.input_nodes + self.output_nodes)*2)
            j_origin = 2*len(self.input_nodes)
            if n_node:
                line_generator[j_origin + 2*(n_node-1)] = 1
                line_generator[j_origin + 2*(n_node)] = -1
                matrix_generator.append(line_generator)

        # Heat flows equalities
        for n_node in range(len(self.input_nodes + self.output_nodes)):
            line_generator = [0]*(len(self.input_nodes + self.output_nodes)*2)
            line_generator[2*(n_node) + 1] = 1
            matrix_generator.append(line_generator)

        matrix = npy.array(matrix_generator)
        return matrix

class NodeEquivalence:
    """
    Defines a link between 3 thermal nodes
    """
    def __init__(self, nodes):
        self.nodes = nodes
        self.system_matrix = self.SystemMatrix()
        self.n_equations, self.n_variables = self.system_matrix.shape

    def SystemMatrix(self):
        """
        Computes block equations system
        T1 = T2
        Phi1 + Phi2 = 0
        """
        matrix = npy.array([[1, 0, -1, 0]])
        return matrix

class ThermalPipe:
    """
    Defines a pipe
    """
    def __init__(self, nodes):
        self.nodes = nodes
        self.system_matrix = self.SystemMatrix()
        self.n_equations, self.n_variables = self.system_matrix.shape

    def SystemMatrix(self, flow=[1], fluid=water):
        """
        TODO Docstring
        """
        constant = fluid.rho*fluid.heat_capacity*flow[0]
        matrix = npy.array([[0.5, 0, 0.5, 0, -1, 0],
                            [0, 1, 0, 0, 0, 0],
                            [0, 0, 0, 1, 0, 0],
                            [constant, 0, -constant, 0, 0, 1]])

        return matrix

class TemperatureBound:
    """
    Defines bounds conditions
    """
    def __init__(self, nodes, value):
        self.nodes = nodes
        self.value = value
        self.system_matrix = self.SystemMatrix()
        self.n_equations, self.n_variables = self.system_matrix.shape

    def SystemMatrix(self):
        """
        TODO Docstring
        """
        system_matrix = npy.array([[1, 0]])
        return system_matrix

class HeatFlowInBound:
    """
    Defines bounds conditions
    """
    def __init__(self, nodes, value):
        self.nodes = nodes
        self.value = value
        self.system_matrix = self.SystemMatrix()
        self.n_equations, self.n_variables = self.system_matrix.shape

    def SystemMatrix(self):
        """
        TODO Docstring
        """
        system_matrix = npy.array([[0, 1]])
        return system_matrix

class HeatFlowOutBound:
    """
    Defines bounds conditions
    """
    def __init__(self, nodes):
        self.nodes = nodes
        self.system_matrix = self.SystemMatrix()
        self.n_equations = 0

    def SystemMatrix(self):
        """
        TODO Docstring
        """
        return []

class Circuit:
    """
    Defines a thermal circuit
    """
    def __init__(self, nodes, blocks):
        self.nodes = nodes
        self.blocks = blocks
        self.temperature_dict, self.flows_dict, self.equations_dict = self.NumberElements()
        self.nodes2blocks = {node : [] for node in nodes}
        for block in blocks:
            for node in block.nodes:
                if block not in self.nodes2blocks[node]:
                    self.nodes2blocks[node].append(block)

    def NumberElements(self):
        temperatures_dict = {}
        flows_dict = {}
        equations_dict = {}
        couples = []
        n_lines = 0
        for j, node in enumerate(self.nodes):
            temperatures_dict[node] = j
            for block in self.blocks:
                if node in block.nodes and (block, node) not in couples:
                    couples.append((block, node))
        j_origin = j + 1
        for j, couple in enumerate(couples):
            j_real = j_origin + j
            flows_dict[couple] = j_real

        for block in self.blocks:
            equations_dict[block] = range(n_lines, n_lines + block.n_equations)
            n_lines += block.n_equations

        return temperatures_dict, flows_dict, equations_dict

    def SystemMatrix(self, input_flows={}, fluid=None):
        """
        Loops on every blocks to get its equation and then assemble the whole system
        """
        n_equations = sum([elt.n_equations for elt in self.blocks + self.nodes])
        n_variables = len(self.temperature_dict.keys()) + len(self.flows_dict.keys())
        system_matrix = npy.zeros((n_equations, n_variables))

        for block in self.blocks:
            if block in input_flows.keys() and len(input_flows[block]) == 1:
                block_system_matrix = block.SystemMatrix(input_flows[block], fluid)
            elif block in input_flows.keys() and len(input_flows[block]) > 1:
                block_system_matrix = block.SystemMatrix(input_flows[block])
            else:
                block_system_matrix = block.SystemMatrix()
            for i_local, i_global in enumerate(self.equations_dict[block]):
                for j, node in enumerate(block.nodes):
                    j_temp_local = 2*j
                    j_hf_local = 2*j + 1
                    j_temp_global = self.temperature_dict[node]
                    j_hf_global = self.flows_dict[(block, node)]
                    system_matrix[i_global][j_temp_global] = block_system_matrix[i_local][j_temp_local]
                    system_matrix[i_global][j_hf_global] = block_system_matrix[i_local][j_hf_local]

        i = i_global + 1
        for node in self.nodes:
            blocks = self.nodes2blocks[node]
            for block in blocks:
                j_hf_global = self.flows_dict[(block, node)]
                system_matrix[i][j_hf_global] = 1
            i += 1

        return system_matrix

    def Solve(self, input_flows={}, fluid=None):
        system_matrix = self.SystemMatrix(input_flows, fluid)
        vector_b = npy.zeros(system_matrix.shape[0])
        bounds_blocks = [block for block in self.blocks\
                         if block.__class__ == TemperatureBound\
                         or block.__class__ == HeatFlowInBound]
        for block in bounds_blocks:
            i_global = self.equations_dict[block][0]
            vector_b[i_global] = block.value
        solution = solve(system_matrix, vector_b)

        return system_matrix, vector_b, solution

    def Draw(self):
        graph = nx.Graph()
        if self.nodes:
            graph.add_nodes_from(self.nodes)
            for block in self.blocks:
                if len(block.nodes) > 1:
                    graph.add_edge(*block.nodes)

        nx.draw_kamada_kawai(graph)

def EquationsSystemAnalysis(M, targeted_variables,
                            arret_surcontrainte=True,
                            dependences=None):
    """
    Analyze an Equation system given by its occurence matrix and the targeted variable to solve
    :return: False (système non solvable si il existe des parties surcontraintes)
    si arret_surcontrainte==True
    sinon, renvoie True, les variables résolvables et l'ordre de résolution
    """
    neq, nvar = M.shape
    G = nx.Graph()
    Gp = nx.DiGraph()
    pos = {}
    for i in range(nvar):
        G.add_node('v'+str(i), bipartite=0)
        Gp.add_node('v'+str(i), bipartite=0)
        pos['v'+str(i)] = [i, 0]

    for i in range(neq):
        G.add_node('e'+str(i), bipartite=1)
        Gp.add_node('e'+str(i), bipartite=1)
        pos['e'+str(i)] = [i, 1]
        for j in range(nvar):
            if M[i, j] != 0:
                G.add_edge('e'+str(i), 'v'+str(j))
                Gp.add_edge('e'+str(i), 'v'+str(j))

    # Adding dependences
    if dependences is not None:
        for eq, var in dependences:
            Gp.add_edge('e'+str(eq), 'v'+str(var))

    sinks = []
    sources = []

    for Gi in nx.connected_component_subgraphs(G):
        M = nx.bipartite.maximum_matching(Gi)

        for n1, n2 in M.items():
            Gp.add_edge(n1, n2)

    for node in Gp.nodes():
        if Gp.out_degree(node) == 0:
            sinks.append(node)
        elif Gp.in_degree(node) == 0:
            sources.append(node)

    G2 = sources[:]
    for node in sources:
        for node2 in nx.descendants(Gp, node):
            if node2 not in G2:
                G2.append(node2)

    if arret_surcontrainte:
        if G2 != []:
            return (False, [], None)

    G3 = sinks[:]
    for node in sinks:
        for node2 in nx.ancestors(Gp, node):
            if node2 not in G3:
                G3.append(node2)

    solvable_vars = []
    for var in targeted_variables:
        if not 'v{}'.format(var) in G2+G3:
            solvable_vars.append(var)

    G1 = G.copy()
    G1.remove_nodes_from(G2+G3)

    G1p = nx.DiGraph()

    G1p.add_nodes_from(G1.nodes())
    for e in G1.edges():
        # Equation vers variable
        if e[0][0] == 'v':
            G1p.add_edge(e[0], e[1])
        else:
            G1p.add_edge(e[1], e[0])


    for G1i in nx.connected_component_subgraphs(G1):
        M1 = nx.bipartite.maximum_matching(G1i)


        for n1, n2 in M1.items():
            if n1[0] == 'e':
                G1p.add_edge(n1, n2)
            else:
                G1p.add_edge(n2, n1)

    # Adding dependences
    if dependences is not None:
        for eq, var in dependences:
            G1p.add_edge('e'+str(eq), 'v'+str(var))

    scc = list(nx.strongly_connected_components(G1p))

    if scc != []:
        C = nx.condensation(G1p, scc)
        # On cherche les indices des blocs contenant chacune des variables
        isc_vars = []
        for isc, sc in enumerate(scc):
            for var in solvable_vars:
                if 'v'+str(var) in sc:
                    isc_vars.append(isc)
                    break
        ancetres_vars = isc_vars[:]

        for isc_var in isc_vars:
            for ancetre in nx.ancestors(C, isc_var):
                if ancetre not in ancetres_vars:
                    ancetres_vars.append(ancetre)

        ordre_sc = [sc for sc in nx.topological_sort(C) if sc in ancetres_vars]
        ordre_ev = []
        for isc in ordre_sc:
            evs = sorted(scc[isc]) # liste d'équations et de variables triées pour être séparées

            levs = int(len(evs)/2)
            ordre_ev.append(([int(e[1:]) for e in evs[0:levs]], [int(v[1:]) for v in evs[levs:]]))

        return (True, solvable_vars, ordre_ev)

    return (False, [], None)

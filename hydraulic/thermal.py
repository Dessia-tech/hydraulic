#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:54:15 2018

@author: jezequel
"""
import numpy as npy
import networkx as nx
from scipy.linalg import solve
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from hydraulic.fluids import water

class Node:
    """
    Defines a thermal node
    """
    def __init__(self, name=''):
        self.n_equations = 1
        self.name = name

    def __repr__(self):
        return self.name

class Resistor:
    """
    Defines a thermal resistor
    """
    def __init__(self, nodes, thickness, thermal_conductivity, area, name=''):
        self.nodes = nodes
        self.thickness = thickness
        self.thermal_conductivity = thermal_conductivity
        self.area = area
        self.resistance = thickness/(thermal_conductivity*area)
        self.system_matrix = self.SystemMatrix()
        self.n_equations, self.n_variables = self.system_matrix.shape
        self.name = name

    def __repr__(self):
        return self.name

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
    def __init__(self, input_nodes, output_nodes, name=''):
        self.input_nodes = input_nodes
        self.output_nodes = output_nodes
        self.nodes = input_nodes + output_nodes
        input_flows = [1 for i in range(len(input_nodes))]
        self.system_matrix = self.SystemMatrix(input_flows)
        self.n_equations, self.n_variables = self.system_matrix.shape
        self.name = name

    def __repr__(self):
        str0 = str(len(self.input_nodes))
        str1 = str(len(self.output_nodes))
        return str0+"-junct-"+str1

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

    def UpdateNode(self, old_node, new_node):
        """
        Replaces old_node by new_node
        """
        index_global = self.nodes.index(old_node)
        if index_global < len(self.input_nodes):
            index = self.input_nodes.index(old_node)
            self.input_nodes[index] = new_node
        else:
            index = self.output_nodes.index(old_node)
            self.output_nodes[index] = new_node
        self.nodes[index_global] = new_node

class NodeEquivalence:
    """
    Defines a link between 3 thermal nodes
    """
    def __init__(self, nodes):
        self.nodes = nodes
        self.system_matrix = self.SystemMatrix()
        self.n_equations, self.n_variables = self.system_matrix.shape

        str0 = self.nodes[0].name
        str1 = self.nodes[1].name
        self.name = str0+"-nodeq-"+str1

    def __repr__(self):
        return self.name

    def SystemMatrix(self):
        """
        Computes block equations system
        T1 = T2
        """
        matrix = npy.array([[1, 0, -1, 0]])
        return matrix

class ThermalPipe:
    """
    Defines a pipe
    """
    def __init__(self, nodes, name=''):
        self.nodes = nodes
        self.system_matrix = self.SystemMatrix()
        self.n_equations, self.n_variables = self.system_matrix.shape
        self.name = name

    def __repr__(self):
        str0 = self.nodes[0].name
        str1 = self.nodes[1].name
        return str0+"-tpipe-"+str1

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

    def UpdateNode(self, old_node, new_node):
        """
        Replaces old_node by new_node
        """
        index = self.nodes.index(old_node)
        self.nodes[index] = new_node

class TemperatureBound:
    """
    Defines bounds conditions
    """
    def __init__(self, nodes, value, name=''):
        self.nodes = nodes
        self.value = value
        self.system_matrix = self.SystemMatrix()
        self.n_equations, self.n_variables = self.system_matrix.shape
        self.name = name

    def __repr__(self):
        return self.name

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
    def __init__(self, nodes, value, name=''):
        self.nodes = nodes
        self.value = value
        self.system_matrix = self.SystemMatrix()
        self.n_equations, self.n_variables = self.system_matrix.shape
        self.name = name

    def __repr__(self):
        return self.name

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
    def __init__(self, nodes, name=''):
        self.nodes = nodes
        self.system_matrix = self.SystemMatrix()
        self.n_equations = 0
        self.name = name

    def __repr__(self):
        return self.name

    def SystemMatrix(self):
        """
        TODO Docstring
        """
        return []

class Circuit:
    """
    Defines a thermal circuit
    """
    def __init__(self, nodes, blocks, name=''):
        self.nodes = nodes
        self.blocks = blocks
        self.nodes2blocks = {node : [] for node in nodes}
        for block in blocks:
            for node in block.nodes:
                if block not in self.nodes2blocks[node]:
                    self.nodes2blocks[node].append(block)
        self.AddEquivalenceBlocks()

        self.temperature_dict, self.flows_dict, self.equations_dict = self.NumberElements()
        self.graph = self.GenerateGraph()
        self.name = name

    def NumberElements(self):
        """
        Gives each variable its global position in solution vector and system matrix
        """
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

    def AddEquivalenceBlocks(self):
        """
        Finds nodes that are connected to 2 hydraulic blocks and replace them
        with node equivalences
        """
        for node in self.nodes:
            hydraulic_blocks = self.nodes2blocks[node]
            are_hydraulic = [IsHydraulicBlock(block) for block in hydraulic_blocks]
            if all(are_hydraulic) and len(hydraulic_blocks) > 1:
                block1 = hydraulic_blocks[0]
                block2 = hydraulic_blocks[1]
                output_node = Node('node'+str(len(self.nodes)))
                equivalence_block = NodeEquivalence([node, output_node])
                block2.UpdateNode(node, output_node)
                self.nodes.append(output_node)
                self.blocks.append(equivalence_block)
                self.nodes2blocks[node] = [block1, equivalence_block]
                self.nodes2blocks[output_node] = [equivalence_block, block2]

    def SystemMatrix(self, input_flows={}, fluid=None):
        """
        Loops on every blocks to get its equation and then assemble the whole system
        """
        n_equations = sum([elt.n_equations for elt in self.blocks + self.nodes])
        n_variables = len(self.temperature_dict.keys()) + len(self.flows_dict.keys())
        system_matrix = npy.zeros((n_equations, n_variables))

        for block in self.blocks:
            if block.__class__ == ThermalPipe:
                block_system_matrix = block.SystemMatrix(input_flows[block], fluid)
            elif block.__class__ == Junction:
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
        """
        Solves system matrix with given input flows and fluid
        """
        system_matrix = self.SystemMatrix(input_flows, fluid)
        vector_b = npy.zeros(system_matrix.shape[0])
        bounds_blocks = [block for block in self.blocks\
                         if block.__class__ == TemperatureBound\
                         or block.__class__ == HeatFlowInBound]
        for block in bounds_blocks:
            i_global = self.equations_dict[block][0]
            vector_b[i_global] = block.value
        solution = solve(system_matrix, vector_b)

        result = ThermalResult(self, system_matrix, vector_b, solution)
        return result

    def GenerateGraph(self):
        """
        Generate circuit graph
        """
        graph = nx.Graph()
        if self.nodes:
            graph.add_nodes_from(self.nodes, node_type='node', node_shape='o')

        if self.blocks:
            graph.add_nodes_from(self.blocks, node_type='block', node_shape='s')

        for block in self.blocks:
            for node in block.nodes:
                graph.add_edge(block, node)

        return graph

    def Draw(self):
        """
        Draws circuit graph with kamada-kawai layout
        """
        nx.draw_kamada_kawai(self.graph)

class ThermalResult:
    """
    Defines a thermal circuit solution
    """
    def __init__(self, circuit, system_matrix, vector_b, solution):
        self.circuit = circuit
        self.system_matrix = system_matrix
        self.vector_b = vector_b
        self.solution = solution

    def Display(self, color_map='jet'):
        """
        Displays solution as a kamada-kawai graph layout
        """
        graph = self.circuit.graph

        graph_pos = nx.kamada_kawai_layout(graph)
        fig = plt.figure("Solution")
        ax = fig.add_subplot(111)

        nodes = [node for node in graph.nodes if graph.nodes[node]['node_type'] == 'node']
        blocks = [node for node in graph.nodes if graph.nodes[node]['node_type'] == 'block']
        block_node_couples = [(block, node) for block in blocks for node in block.nodes]

        temp_values = [self.solution[self.circuit.temperature_dict[node]] for node in nodes]

        temp_pts = nx.draw_networkx_nodes(graph,
                                          pos=graph_pos,
                                          nodelist=nodes,
                                          node_color=temp_values,
                                          node_size=100,
                                          cmap=color_map)
        nx.draw_networkx_nodes(graph,
                               pos=graph_pos,
                               nodelist=blocks,
                               node_shape='s',
                               node_color='k',
                               node_size=50)

        for bnc in block_node_couples:
            index = self.circuit.flows_dict[bnc]
            flow_value = self.solution[index]
            if flow_value <= 0:
                # From block to node
                edge_position = (graph_pos[bnc[0]], graph_pos[bnc[1]])
            else:
                # From node to block
                edge_position = (graph_pos[bnc[1]], graph_pos[bnc[0]])

            x = edge_position[0][0]
            y = edge_position[0][1]
            dx = edge_position[1][0] - x
            dy = edge_position[1][1] - y
            ax.add_patch(patches.FancyArrow(x, y, dx, dy, width=abs(flow_value)/1000,
                                            length_includes_head=True,
                                            head_width=abs(flow_value)*1.5/1000))
        ax.axis('equal')
        plt.colorbar(temp_pts)

def IsHydraulicBlock(block):
    """
    Checks block type. Return False if block is only thermal, True otherwise
    """
    if block.__class__ == ThermalPipe or block.__class__ == Junction:
        return True
    return False

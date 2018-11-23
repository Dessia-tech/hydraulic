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
#        self.n_equations = 1
        self.name = name

#    def __repr__(self):
#        return self.name

class Block:
    """
    Defines a thermal block
    """
    def __init__(self, nodes, active_nodes, name=''):
        self.nodes = nodes
        self.active_nodes = active_nodes
        self.name = name

    def __repr__(self):
        return self.name

    def _get_n_equations(self):
        n_equations = self.system_matrix.shape[0]
        return n_equations
    n_equations = property(_get_n_equations)

class Resistor(Block):
    """
    Defines a thermal resistor
    """
    def __init__(self, nodes, thickness, thermal_conductivity, area, name=''):
        self.thickness = thickness
        self.thermal_conductivity = thermal_conductivity
        self.area = area
        self.resistance = thickness/(thermal_conductivity*area)

        self.system_matrix = self.SystemMatrix()
        Block.__init__(self, nodes, nodes, name)

    def SystemMatrix(self):
        """
        Computes node equations system
        sum(phi) = 0
        """
        matrix = npy.array([[-self.resistance, 1, self.resistance, 0],
                            [0, 1, 0, 1]])
        return matrix

class ThermalPipe(Block):
    """
    Defines a pipe
    """
    def __init__(self, nodes, flow, fluid=water, name=''):
        constant = fluid.rho*fluid.heat_capacity*flow
        matrix_generator = [[0.5, 0, 0.5, 0, -1, 0],
                            [constant, 0, -constant, 0, 0, 1]]
#        if nodes[0] in active_nodes:
#            matrix_generator.append([0, 1, 0, 0, 0, 0])
#        if nodes[1] in active_nodes:
#            matrix_generator.append([0, 0, 0, 1, 0, 0])
        self.system_matrix = npy.array(matrix_generator)
        self.flow = flow
        self.fluid = fluid
        Block.__init__(self, nodes, [nodes[2]], name)

    def __repr__(self):
        str0 = self.nodes[0].name
        str1 = self.nodes[1].name
        return str0+"-tpipe-"+str1

    def SystemMatrix(self):
        """
        TODO Docstring
        """
        constant = self.fluid.rho*self.fluid.heat_capacity*self.flow
        matrix_generator = [[0.5, 0, 0.5, 0, -1, 0],
                            [constant, 0, -constant, 0, 0, 1]]
#        if self.nodes[0] in self.active_nodes:
#            matrix_generator.append([0, 1, 0, 0, 0, 0])
#        if self.nodes[1] in self.active_nodes:
#            matrix_generator.append([0, 0, 0, 1, 0, 0])
        matrix = npy.array(matrix_generator)
        return matrix

    def UpdateFlow(self, new_flow):
        self.flow = new_flow
        self.system_matrix = self.SystemMatrix()

    def UpdateFluid(self, new_fluid):
        self.fluid = new_fluid
        self.system_matrix = self.SystemMatrix()

class Junction(Block):
    """
    Defines a junction, with several inputs and outputs
    """
    def __init__(self, input_nodes, output_nodes,
                 input_flows, fluid=water,
                 name=''):
        self.input_nodes = input_nodes
        self.output_nodes = output_nodes

        self.input_flows = input_flows
        self.fluid = fluid
        self.system_matrix = self.SystemMatrix()
        nodes = input_nodes + output_nodes
        Block.__init__(self, nodes, [], name)

    def __repr__(self):
        str0 = str(len(self.input_nodes))
        str1 = str(len(self.output_nodes))
        return str0+"-junct-"+str1

    def SystemMatrix(self):
        """
        Computes block equations system
        Output temperatures are a weighted mean of input volumetric flows
        Heat flows are equals
        """
        matrix_generator = []
        line_generator = []

        # Line 1 : Weigthed mean volumetric flows equation
        for n_node in range(len(self.input_nodes)):
            flow = self.input_flows[n_node]
            line_generator.extend([flow/sum(self.input_flows), 0])

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

#        # Heat flows equalities
#        for n_node in range(len(self.input_nodes + self.output_nodes)):
#            line_generator = [0]*(len(self.input_nodes + self.output_nodes)*2)
#            line_generator[2*(n_node) + 1] = 1
#            matrix_generator.append(line_generator)

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

    def UpdateFlow(self, new_flow):
        self.flow = new_flow
        self.system_matrix = self.SystemMatrix()

    def UpdateFluid(self, new_fluid):
        self.fluid = new_fluid
        self.system_matrix = self.SystemMatrix()

class NodeEquivalence(Block):
    """
    Defines a link between 3 thermal nodes
    """
    def __init__(self, nodes):
        self.system_matrix = self.SystemMatrix()
        str0 = nodes[0].name
        str1 = nodes[1].name
        name = str0+"-nodeq-"+str1
        Block.__init__(self, nodes, nodes, name)

    def SystemMatrix(self):
        """
        Computes block equations system
        T1 = T2
        """
        matrix = npy.array([[1, 0, -1, 0]])
        return matrix

    def UpdateNode(self, old_node, new_node):
        """
        Replaces old_node by new_node
        """
        index = self.nodes.index(old_node)
        self.nodes[index] = new_node

class TemperatureBound(Block):
    """
    Defines bounds conditions
    """
    def __init__(self, nodes, value, name=''):
        self.value = value

        self.system_matrix = self.SystemMatrix()
        Block.__init__(self, nodes, nodes, name)

    def SystemMatrix(self):
        """
        TODO Docstring
        """
        system_matrix = npy.array([[1, 0]])
        return system_matrix

class HeatFlowInBound(Block):
    """
    Defines bounds conditions
    """
    def __init__(self, nodes, value, name=''):
        self.value = value

        self.system_matrix = self.SystemMatrix()
        Block.__init__(self, nodes, nodes, name)

    def SystemMatrix(self):
        """
        TODO Docstring
        """
        system_matrix = npy.array([[0, 1]])
        return system_matrix

class HeatFlowOutBound(Block):
    """
    Defines bounds conditions
    """
    def __init__(self, nodes, name=''):
        self.system_matrix = self.SystemMatrix()
        Block.__init__(self, nodes, nodes, name)

    def SystemMatrix(self):
        """
        TODO Docstring
        """
        return npy.empty((0, 0))

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
#        self.AddEquivalenceBlocks()

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
                if node in block.active_nodes and (block, node) not in couples:
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

    def SystemMatrix(self):
        """
        Loops on every blocks to get its equation and then assemble the whole system
        """
        # Initialize system matrix
        active_nodes = []
        for block in self.blocks:
            for node in block.active_nodes:
                if node not in active_nodes:
                    active_nodes.append(node)
        n_equations = sum([block.n_equations for block in self.blocks]) + len(active_nodes)
        n_variables = len(self.temperature_dict.keys()) + len(self.flows_dict.keys())
        system_matrix = npy.zeros((n_equations, n_variables))

        # Write blocks equations
        for block in self.blocks:
            block_system_matrix = block.SystemMatrix()
            for i_local, i_global in enumerate(self.equations_dict[block]):
                for j, node in enumerate(block.nodes):
                    j_temp_local = 2*j
                    j_temp_global = self.temperature_dict[node]
                    system_matrix[i_global][j_temp_global] = block_system_matrix[i_local][j_temp_local]
                    if node in block.active_nodes:
                        j_hf_local = 2*j + 1
                        j_hf_global = self.flows_dict[(block, node)]
                        system_matrix[i_global][j_hf_global] = block_system_matrix[i_local][j_hf_local]

        # Write nodes equations
        i = i_global + 1
        valid = False
        for node in self.nodes:
            blocks = self.nodes2blocks[node]
            for block in blocks:
                if node in block.active_nodes:
                    j_hf_global = self.flows_dict[(block, node)]
                    system_matrix[i][j_hf_global] = 1
                    valid = True
            if valid:
                i += 1
            valid = False

        return system_matrix

    def Solve(self):
        """
        Solves system matrix with given input flows and fluid
        """
        #  Build system matrix
        system_matrix = self.SystemMatrix()

        # Check lines and columns
        vector_b = npy.zeros(system_matrix.shape[0])

        # Write boundary conditions
        bounds_blocks = [block for block in self.blocks\
                         if block.__class__ == TemperatureBound\
                         or block.__class__ == HeatFlowInBound]
        for block in bounds_blocks:
            i_global = self.equations_dict[block][0]
            vector_b[i_global] = block.value

        # Solve
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
        fig = plt.figure()
        ax = fig.add_subplot(111)

        nodes = [node for node in graph.nodes if graph.nodes[node]['node_type'] == 'node']
        blocks = [node for node in graph.nodes if graph.nodes[node]['node_type'] == 'block']
        block_node_couples = [(block, node) for block in blocks for node in block.nodes]

        for bnc in block_node_couples:
            add_link = False
            if bnc in self.circuit.flows_dict.keys():
                index = self.circuit.flows_dict[bnc]
                flow_value = self.solution[index]
                if flow_value < 0:
                    # From block to node
                    edge_position = (graph_pos[bnc[0]], graph_pos[bnc[1]])
                elif flow_value > 0:
                    # From node to block
                    edge_position = (graph_pos[bnc[1]], graph_pos[bnc[0]])
                else:
                    add_link = True
    
                if not add_link:
                    x = edge_position[0][0]
                    y = edge_position[0][1]
                    dx = edge_position[1][0] - x
                    dy = edge_position[1][1] - y
                    ax.add_patch(patches.FancyArrow(x, y, dx, dy, width=abs(flow_value)/1000,
                                                    length_includes_head=True,
                                                    head_width=abs(flow_value)*1.5/1000))
            else:
                add_link = True

            if add_link:
                edge_position = (graph_pos[bnc[0]], graph_pos[bnc[1]])
                x = [edge_position[0][0], edge_position[1][0]]
                y = [edge_position[0][1], edge_position[1][1]]
                ax.plot(x, y, 'k-')
                
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
        
        ax.axis('equal')
        plt.colorbar(temp_pts)

def IsHydraulicBlock(block):
    """
    Checks block type. Return False if block is only thermal, True otherwise
    """
    if block.__class__ == ThermalPipe or block.__class__ == Junction:
        return True
    return False
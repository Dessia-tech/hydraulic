#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import numpy as npy
import networkx as nx
from scipy.linalg import solve
import matplotlib.pyplot as plt
from matplotlib import patches
from hydraulic.fluids import water
import volmdlr as vm


class Node:
    """
    Defines a thermal node
    """
    def __init__(self, name=''):
        self.name = name

    def __str__(self):
        return self.name


class Block:
    """
    Defines a thermal block. This is an abstract class
    """
    def __init__(self, nodes, active_nodes, name=''):
        self.nodes = nodes
        self.active_nodes = active_nodes
        self.name = name
        
        self._utd_system_equations = False

    def __str__(self):
        return self.name
        
    def _get_system_matrix(self):
        if not self._utd_system_equations:
            self._system_matrix, self._system_rhs = self.SystemEquations()
            self._utd_system_equations = True
        return self._system_matrix
    
    system_matrix = property(_get_system_matrix)

    def _get_system_rhs(self):
        if not self._utd_system_equations:
            self._system_matrix, self._system_rhs = self.SystemEquations()
            self._utd_system_equations = True
        return self._system_rhs

    system_rhs = property(_get_system_rhs)

    def _get_n_equations(self):
        n_equations = self.system_matrix.shape[0]
        return n_equations
    n_equations = property(_get_n_equations)


class Resistor(Block):
    """
    Defines a thermal resistor
    """
    def __init__(self, nodes, resistance, name=''):
        self.resistance = resistance

        Block.__init__(self, nodes, nodes, name)

    def __str__(self):
        str0 = self.nodes[0].name
        str1 = self.nodes[1].name
        return str0+"-res-"+str1

    def SystemEquations(self):
        """
        Computes node equations system
        sum(phi) = 0
        """
        matrix = npy.array([[-1, self.resistance, 1, 0],
                            [0, 1, 0, 1]])
        b = npy.array([0, 0])
        return matrix, b


class ThermalPipe(Block):
    """
    Defines a pipe
    """
    def __init__(self, nodes, flow, fluid=water, name=''):
        self.flow = flow
        self.fluid = fluid
        Block.__init__(self, nodes, [nodes[2]], name)

    def __str__(self):
        str0 = self.nodes[0].name
        str1 = self.nodes[1].name
        return str0+"-tpipe-"+str1

    def SystemEquations(self):
        """
        TODO Docstring
        """
        constant = self.fluid.rho*self.fluid.heat_capacity*self.flow
        matrix = npy.array([[0.5, 0, 0.5, 0, -1, 0],
                            [constant, 0, -constant, 0, 0, 1]])
        b = npy.array([0, 0])
        return matrix, b

    def UpdateFlow(self, new_flow):
        self.flow = new_flow
        self.system_matrix, self.system_rhs = self.SystemEquations()

    def UpdateFluid(self, new_fluid):
        self.fluid = new_fluid
        self.system_matrix, self.system_rhs = self.SystemEquations()


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
        nodes = input_nodes + output_nodes
        Block.__init__(self, nodes, [], name)

    def __str__(self):
        str0 = str(len(self.input_nodes))
        str1 = str(len(self.output_nodes))
        return str0+"-junct-"+str1

    def SystemEquations(self):
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
        lon = len(self.output_nodes)
        for n_node in range(lon):
            line_generator = [0]*(len(self.input_nodes + self.output_nodes)*2)
            j_origin = 2*len(self.input_nodes)
            if n_node:
                line_generator[j_origin + 2*(n_node-1)] = 1
                line_generator[j_origin + 2*(n_node)] = -1
                matrix_generator.append(line_generator)

        matrix = npy.array(matrix_generator)
        b = npy.zeros(lon+1)
        return matrix, b

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
        self.system_matrix = self.SystemEquations()

    def UpdateFluid(self, new_fluid):
        self.fluid = new_fluid
        self.system_matrix = self.SystemEquations()


class NodeEquivalence(Block):
    """
    Defines a link between 3 thermal nodes
    """
    def __init__(self, nodes):
        str0 = nodes[0].name
        str1 = nodes[1].name
        name = str0+"-nodeq-"+str1
        Block.__init__(self, nodes, nodes, name)

    def SystemEquations(self):
        """
        Computes block equations system
        T1 = T2
        """
        matrix = npy.array([[1, 0, -1, 0],
                            [0, 1, 0, 1]])
        b = npy.array([0, 0])
        return matrix, b

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

        Block.__init__(self, nodes, nodes, name)

    def SystemEquations(self):
        """
        TODO Docstring
        """
        system_matrix = npy.array([[1, 0]])
        rhs = npy.array([self.value])
        return system_matrix, rhs


class HeatFlowInBound(Block):
    """
    Defines bounds conditions
    """
    def __init__(self, nodes, value, name=''):
        self.value = value

        Block.__init__(self, nodes, nodes, name)

    def SystemEquations(self):
        """
        TODO Docstring
        """
        system_matrix = npy.array([[0, 1]])
        rhs = npy.array([self.value])
        return system_matrix, rhs


class HeatFlowOutBound(Block):
    """
    Defines bounds conditions
    """
    def __init__(self, nodes, name=''):
        Block.__init__(self, nodes, nodes, name)

    def SystemEquations(self):
        """
        TODO Docstring
        """
        return npy.empty((0, 0)), npy.empty(0)
    
class UnidimensionalMedium(Block):
    """
    Defines an unidimensional medium of length l
    """
    def __init__(self, nodes, conductivity, length, contact_area, name=''):
#        self.resistance = resistance_factor/area_factor

        self.conductivity = conductivity
        self.length = length
        self.contact_area = contact_area
        Block.__init__(self, nodes, nodes, name)

    def __str__(self):
        str0 = self.nodes[0].name
        str1 = self.nodes[1].name
        return str0+"-res-"+str1

    def SystemEquations(self):
        """
        Computes node equations system
        sum(phi) = 0
        """
        matrix = npy.array([[0, 1, 0, 1, 0, 1],
                            [-1, self.length/self.conductivity/self.contact_area, 1, 0, 0, self.length*0.5/self.conductivity/self.contact_area],
                            [-0.5, 0, -0.5, 0, 1, 0]])
        rhs = npy.zeros(3)
    
        return matrix, rhs


class Circuit:
    """
    Defines a thermal circuit
    """
    def __init__(self, nodes, blocks, name=''):
        self.nodes = nodes
        self.blocks = blocks
        self.nodes2blocks = self.SetNodes2Blocks()

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

#    def AddEquivalenceBlocks(self):
#        """
#        Finds nodes that are connected to 2 hydraulic blocks and replace them
#        with node equivalences
#        """
#        for node in self.nodes:
#            hydraulic_blocks = self.nodes2blocks[node]
#            are_hydraulic = [IsHydraulicBlock(block) for block in hydraulic_blocks]
#            if all(are_hydraulic) and len(hydraulic_blocks) > 1:
#                block1 = hydraulic_blocks[0]
#                block2 = hydraulic_blocks[1]
#                output_node = Node('node'+str(len(self.nodes)))
#                equivalence_block = NodeEquivalence([node, output_node])
#                block2.UpdateNode(node, output_node)
#                self.nodes.append(output_node)
#                self.blocks.append(equivalence_block)
#                self.Update()
##                self.nodes2blocks[node] = [block1, equivalence_block]
##                self.nodes2blocks[output_node] = [equivalence_block, block2]

    def AddBlocks(self, blocks):
        for block in blocks:
            self.blocks.append(block)
            for node in block.nodes:
                if node not in self.nodes:
                    self.nodes.append(node)
        self.Update()

    def AddResistor(self, node,
                    resistance_factor, area_factor,
                    boundary_condition):
        bc_node = boundary_condition.nodes[0]
        resistor = Resistor([node, bc_node],
                            resistance_factor, area_factor)
        self.nodes.append(bc_node)
        self.blocks.extend([resistor, boundary_condition])
        self.Update()

    def Update(self):
        self.nodes2blocks = self.SetNodes2Blocks()
        self.temperature_dict, self.flows_dict, self.equations_dict = self.NumberElements()
        self.graph = self.GenerateGraph()

    def SetNodes2Blocks(self):
        nodes2blocks = {node : [] for node in self.nodes}
        for block in self.blocks:
            for node in block.nodes:
                if block not in nodes2blocks[node]:
                    nodes2blocks[node].append(block)
        return nodes2blocks

    def SystemEquations(self):
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
        rhs = npy.zeros(n_equations)

        # Write blocks equations
        for block in self.blocks:
            matrix_block, rhs_block = block.SystemEquations()
            for i_local, i_global in enumerate(self.equations_dict[block]):
                for j, node in enumerate(block.nodes):
                    j_temp_local = 2*j
                    j_temp_global = self.temperature_dict[node]
                    system_matrix[i_global][j_temp_global] = matrix_block[i_local][j_temp_local]
                    rhs[i_global] = rhs_block[i_local]
                    if node in block.active_nodes:
                        j_hf_local = 2*j + 1
                        j_hf_global = self.flows_dict[(block, node)]
                        system_matrix[i_global][j_hf_global] = matrix_block[i_local][j_hf_local]
                        rhs[i_global] = rhs_block[i_local]
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

        return system_matrix, rhs

    def Solve(self):
        """
        Solves system matrix with given input flows and fluid
        """
        #  Build system matrix
        system_matrix, vector_b = self.SystemEquations()

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
        plt.figure()
        graph_pos = nx.kamada_kawai_layout(self.graph)

        nodes = [node for node in self.graph.nodes
                 if self.graph.nodes[node]['node_type'] == 'node']
        blocks = [node for node in self.graph.nodes
                  if self.graph.nodes[node]['node_type'] == 'block']

        nx.draw_networkx_nodes(self.graph,
                               pos=graph_pos,
                               nodelist=nodes,
                               node_color='r',
                               node_size=100)
        nx.draw_networkx_nodes(self.graph,
                               pos=graph_pos,
                               nodelist=blocks,
                               node_shape='s',
                               node_color='k',
                               node_size=50)
        nx.draw_networkx_edges(self.graph, graph_pos)


class ThermalResult:
    """
    Defines a thermal circuit solution
    """
    def __init__(self, circuit, system_matrix, vector_b, solution):
        self.circuit = circuit
        self.system_matrix = system_matrix
        self.vector_b = vector_b
        self.solution = solution

    def Display(self, ax=None, color_map='jet', max_width=0.01,
                position=None, x3D=vm.x3D, y3D=vm.y3D):
        """
        Displays solution as a kamada-kawai graph layout
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        graph = self.circuit.graph

        nodes = [node for node in graph.nodes() if graph.nodes[node]['node_type'] == 'node']
        blocks = [node for node in graph.nodes() if graph.nodes[node]['node_type'] == 'block']

        if position is None:
            graph_pos = nx.kamada_kawai_layout(graph)

            block_node_couples = [(block, node) for block in blocks for node in block.nodes]
            heat_flows = {}
            for bnc in block_node_couples:
                if bnc in self.circuit.flows_dict.keys():
                    index = self.circuit.flows_dict[bnc]
                    heat_flows[bnc] = self.solution[index]
                else:
                    edge_position = (graph_pos[bnc[0]], graph_pos[bnc[1]])
                    x = [edge_position[0][0], edge_position[1][0]]
                    y = [edge_position[0][1], edge_position[1][1]]
                    ax.plot(x, y, 'k-')
    
            width = float(max_width/max(heat_flows.values()))
    
            for bnc, heat_flow in heat_flows.items():
                add_link = False
                if heat_flow < 0:
                    # From block to node
                    edge_position = (graph_pos[bnc[0]], graph_pos[bnc[1]])
                elif heat_flow > 0:
                    # From node to block
                    edge_position = (graph_pos[bnc[1]], graph_pos[bnc[0]])
                else:
                    add_link = True
    
                if add_link:
                    edge_position = (graph_pos[bnc[0]], graph_pos[bnc[1]])
                    x = [edge_position[0][0], edge_position[1][0]]
                    y = [edge_position[0][1], edge_position[1][1]]
                    ax.plot(x, y, 'k-')
                else:
                    x = edge_position[0][0]
                    y = edge_position[0][1]
                    dx = edge_position[1][0] - x
                    dy = edge_position[1][1] - y
                    ax.add_patch(patches.FancyArrow(x, y, dx, dy, width=abs(width*heat_flow),
                                                    length_includes_head=True,
                                                    head_width=abs(width*heat_flow)*2,
                                                    edgecolor='k',
                                                    facecolor='gray'))

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
                                   node_size=10)
        else:
            graph_pos = {node : (point.Dot(x3D), point.Dot(y3D))
                         for node, point in position.items()}
            nodes = graph_pos.keys()
            
            # Blocks positions : centroid of points
            block_pos = {}
            for block in list(blocks):
                if IsHydraulicBlock(block):
                    coords = npy.asarray([graph_pos[n] for n in block.nodes if n in nodes])
                    length = coords.shape[0]
                    sum_x = npy.sum(coords[:, 0])
                    sum_y = npy.sum(coords[:, 1])
                    block_pos[block] = npy.array((sum_x/length, sum_y/length))

            edges = [edge for edge in graph.edges() if edge[0] in nodes and edge[1] in block_pos.keys()]

            temp_values = [self.solution[self.circuit.temperature_dict[node]] for node in nodes]
    
            temp_pts = nx.draw_networkx_nodes(graph,
                                              pos=graph_pos,
                                              nodelist=nodes,
                                              node_color=temp_values,
                                              node_size=100,
                                              cmap=color_map)
            nx.draw_networkx_nodes(graph,
                                   pos=block_pos,
                                   nodelist=block_pos.keys(),
                                   node_shape='s',
                                   node_color='k',
                                   node_size=10)
            
            graph_pos.update(block_pos)
            nx.draw_networkx_edges(graph,
                                   pos=graph_pos,
                                   edgelist=edges)
            

        ax.axis('equal')
        plt.colorbar(temp_pts)


class ThermohydraulicCircuit:
    def __init__(self, hydraulic_circuit, thermal_circuit,
                 interface_nodes,
                 point2node, pipe2block):
        self.hydraulic_circuit = hydraulic_circuit
        self.thermal_circuit = thermal_circuit
        self.interface_nodes = interface_nodes
        self.point2node = point2node
        self.pipe2block = pipe2block

        self.node2point = {}
        self.block2pipe = {}
        for point, node in point2node.items():
            self.node2point[node] = point
        for pipe, block in pipe2block.items():
            self.block2pipe[block] = pipe

    def Draw(self, x3D=vm.x3D, y3D=vm.y3D, position=False):
        """
        TODO Docstring
        """
        if position:
            graph_pos = {node : (point.Dot(x3D), point.Dot(y3D))
                         for node, point in self.node2point.items()}

            thermal_graph = self.thermal_circuit.graph
            nodes = graph_pos.keys()
            
            nx.draw_networkx_nodes(thermal_graph,
                                   pos=graph_pos,
                                   nodelist=nodes,
                                   node_color='r',
                                   node_size=100)

            # Blocks positions : centroid of points
            block_pos = {}
            for pipe in self.hydraulic_circuit.pipes:
                coords = npy.asarray([(point.Dot(x3D), point.Dot(y3D)) for point in pipe.active_points])
                length = coords.shape[0]
                sum_x = npy.sum(coords[:, 0])
                sum_y = npy.sum(coords[:, 1])
                block = self.pipe2block[pipe]
                block_pos[block] = npy.array((sum_x/length, sum_y/length))
            
            blocks = block_pos.keys()
            nx.draw_networkx_nodes(thermal_graph,
                                   pos=block_pos,
                                   nodelist=blocks,
                                   node_shape='s',
                                   node_color='k',
                                   node_size=10)
            
            graph_pos.update(block_pos)
            edges = [edge for edge in thermal_graph.edges if edge[0] in nodes and edge[1] in blocks]
                
            nx.draw_networkx_edges(thermal_graph,
                                   pos=graph_pos,
                                   edgelist=edges)

        else:
            self.thermal_circuit.Draw()


def IsHydraulicBlock(block):
    """
    Checks block type. Return False if block is only thermal, True otherwise
    """
    if block.__class__ == ThermalPipe or block.__class__ == Junction:
        return True
    return False

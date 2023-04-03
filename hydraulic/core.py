#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hydraulic circuit definition
"""

import math
import numpy as npy
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import patches
from scipy.linalg import solve
import hydraulic.pipes as hyp
import hydraulic.thermal as th
from hydraulic import fluids
import volmdlr as vm
import volmdlr.core as vmc
from copy import copy
import dessia_common.core as dcc

# Definition of equivalent L/D values
# Lists in increasing order
ben_angle = [0, 45, 90, 180, 360]
ben_LD = [0, 16, 30, 60, 120]
# Bend radius ?
mbd_angle = [0, 45, 90]
mbd_LD = [0, 15, 58]
enl_rap = [1/4, 1/2, 3/4, 1]
enl_LD = [30, 20, 6.5, 0]
ctr_rap = [1, 4/3, 2, 4]
ctr_LD = [0, 6.5, 11, 15]

class Circuit(dcc.PhysicalObject):
    """
    General class for 2D/3D circuits.
    """
    def __init__(self, points, pipes, boundary_conditions, fluid, name:str=''):
        self.points = points
        self.pipes = pipes
        self.boundary_conditions = boundary_conditions
        self.fluid = fluid
        self.name = name

        self.non_interior_points = {point for pipe in self.pipes\
                                    for point in pipe.active_points}

        self.graph = self.get_graph()
        # Creating thermal nodes out of circuit points

    def get_graph(self):
        graph = nx.Graph()
        if self.points and self.pipes:

            graph.add_nodes_from(self.non_interior_points,
                                 node_type='point',
                                 node_shape='o')
            graph.add_nodes_from(self.pipes,
                                 node_type='pipe',
                                 node_shape='s')
            graph.add_nodes_from(self.boundary_conditions,
                                 node_type='boundary_condition',
                                 node_shape='^')
            
        

        for pipe in self.pipes:
            for point in pipe.active_points:
                graph.add_edge(pipe, point)
                
        for bc in self.boundary_conditions:
            graph.add_edge(bc, bc.points[0])
        
        return graph


    def points_to_pipes(self):
        point2pipes = {point : [] for point in self.points}
        for point in self.points:
            for pipe in self.pipes + self.boundary_conditions:
                if point in pipe.points:
                    point2pipes[point].append(pipe)
        return point2pipes

    def plot_graph(self):
        graph_pos = nx.kamada_kawai_layout(self.graph)
        labels = {node : str(node) for node in list(self.graph.nodes)}
        plt.figure("circuit graph")
        # nx.draw_networkx(self.graph,
        #                  pos=graph_pos,
        #                  labels=labels_dict,
        #                  withlabels=True,
        #                  fontweight="bold")

        nx.draw_networkx_nodes(self.graph,
                               nodelist=self.non_interior_points,
                               pos=graph_pos)
        nx.draw_networkx_nodes(self.graph,
                               nodelist=self.pipes,
                               pos=graph_pos,
                               node_shape='s',
                               node_color='g')
        # Boundary conditions
        nx.draw_networkx_nodes(self.graph,
                               nodelist=self.boundary_conditions,
                               pos=graph_pos,
                               node_shape='h',
                               node_color='y')
        nx.draw_networkx_edges(self.graph,
                               pos=graph_pos)
        nx.draw_networkx_nodes(self.graph,
                               nodelist=self.pipes,
                               pos=graph_pos,
                               node_shape='s',
                               node_color='grey')
        
        nx.draw_networkx_labels(self.graph,
                                pos=graph_pos,
                                labels=labels)
        

    def settings(self):
        pressures_dict = {}
        flows_dict = {}
        equations_dict = {}
        nvars = 0
        neqs = 0

        for point in self.points:
            pressures_dict[point] = nvars
            nvars += 1

        for pipe in self.pipes:
            for point in pipe.active_points:
                flows_dict[(pipe, point)] = nvars
                nvars += 1
                
            equations_dict[pipe] = range(neqs, neqs + pipe.n_equations)
            neqs += pipe.n_equations

        # Counting equations for boundary conditions
        for bc in self.boundary_conditions:
            for point in bc.active_points:
                flows_dict[(bc, point)] = nvars
                nvars += 1
                
            equations_dict[bc] = range(neqs, neqs + bc.n_equations)
            neqs += bc.n_equations
           
        points2pipes = self.points_to_pipes()
        for point in self.points:
            active = False
            pipes = points2pipes[point]
            for pipe in pipes:
                if point in pipe.active_points:
                    active = True
                    break
            if active:
                neqs += 1

        return pressures_dict, flows_dict, equations_dict, neqs, nvars

    def system_matrix(self, constant):
        pressures_dict, flows_dict, equations_dict, neqs, nvars = self.settings()
        system_matrix = npy.zeros((neqs, nvars))
        points2pipes = self.points_to_pipes()

        # Write blocks equations
        for pipe in self.pipes + self.boundary_conditions:
            if pipe in self.boundary_conditions:
                block_system_matrix = pipe.system_matrix()
            else:
                block_system_matrix = pipe.system_matrix(constant)
            for i_local, i_global in enumerate(equations_dict[pipe]):
                for j, point in enumerate(pipe.points):
                    j_press_local = 2*j
                    j_press_global = pressures_dict[point]
                    system_matrix[i_global][j_press_global] = block_system_matrix[i_local][j_press_local]
                    if point in pipe.active_points:
                        j_flow_local = 2*j + 1
                        j_flow_global = flows_dict[(pipe, point)]
                        system_matrix[i_global][j_flow_global] = block_system_matrix[i_local][j_flow_local]

        i = i_global + 1
        valid = False
        for point in self.points:
            pipes = points2pipes[point]
            for pipe in pipes:
                if point in pipe.active_points:
                    j_flow_global = flows_dict[(pipe, point)]
                    system_matrix[i][j_flow_global] = 1
                    valid = True
            if valid:
                i += 1
            valid = False

        return system_matrix

    def solve_fluidics(self):
        """
        solves the hydraulic circuit for
        :param imposed: = "pressure" or "flow"
        pressure : imposed_values = deltaP = dictionary {point:pressure}
        flow : imposed_values = Qs = dictionary {point:flow rate}
               flow is set at input points
               for output points specify "out" as value
               output points will be set at pressure = low_p
        """
        dictionaries = self.settings()
        equations_dict = dictionaries[2]

        constant = 2/(self.fluid.rho*self.fluid.nu)
        system_matrix = self.system_matrix(constant)
        vector_b = npy.zeros(system_matrix.shape[0])

        for boundary_condition in self.boundary_conditions:
            j_global = equations_dict[boundary_condition]
            vector_b[j_global] = boundary_condition.value

        # Solve
        solution = solve(system_matrix, vector_b)
        result = FluidicsResults(self, solution)
        return result

    def VerifyQ(self):
        """
        Check if the flow rate is laminar inside the circuit
        """
        turbulent = 0
        for pipe in self.pipes:
            Re = 2*pipe.flow/(self.fluid.nu*math.pi*pipe.radius)
            if Re > 2000.:  # Laminar condition
                turbulent = 1
                print("Warning, flow is not laminar in {} with Re = {}".format(pipe, Re))
        if not turbulent:
            print("Flow is laminar in the circuit")

    def to_dict(self):
        d = {'points': [point.to_dict() for point in self.points]}
        d['pipes'] = [pipe.to_dict() for pipe in self.pipes]
        d['boundary_conditions'] = [bc.to_dict() for bc in self.boundary_conditions]
        d['fluid'] = self.fluid.to_dict()
        return d


class Circuit2D(Circuit):
    def __init__(self, points, pipes, boundary_conditions, fluid, name: str = ''):
        Circuit.__init__(self, points, pipes, boundary_conditions, fluid, name=name)

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        for pipe in self.pipes:
            pipe.plot(ax)

    def to_gmsh(self, name="Generated_Circuit", path="", size=1., index=0):
        """
        Exports the circuit as a 1D geo file for gmsh
        Size defines the size of the mesh
        """
        file_name = path+name+".geo"
        file = open(file_name, "w")
        n_pts = len(self.points)
        pts_index = dict()
        geo_txt = ""
        for i in range(n_pts):
            point = self.points[i]
            geo_txt += "Point({0}) = [{1[0]},{1[1]},0.,{2}];\n".format(i+1, point, size)
            pts_index[point] = i+1
        n_pipes = len(self.pipes)
        added_elems = 0
        physicals = []  # Index of lines for physical lines
        for j in range(n_pipes):
            pipe = self.pipes[j]
            [mess, phys] = pipe.gmsh_repr1d(j+1+added_elems, pts_index)
            physicals.append(phys)
            geo_txt += mess
            added_elems += len(phys)-1
        geo_txt = geo_txt.replace("[", "{")
        geo_txt = geo_txt.replace("]", "}")
        file.write(geo_txt)
        file.close()
        if index:
            return (pts_index, n_pipes+added_elems, physicals)

    def FindLimits(self):
        # Returns the nodes which are defined as the end of a pipe
        points_lim = []
        for point in self.points:
            if len(self.graph.adj[point]) == 1:
                points_lim.append(point)
        return points_lim


class Circuit3D(Circuit):
    """
    Defines a 3D circuit.

    """

    def __init__(self, points, pipes, boundary_conditions, fluid, name: str = ''):
        Circuit.__init__(self, points, pipes, boundary_conditions, fluid, name=name)

    def plot2d(self, x3D=vm.X3D, y3D=vm.Y3D, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        for pipe in self.pipes:
            pipe.plot2d(x3D, y3D, ax=ax)
        return ax

    def plot(self):
        ax = self.pipes[0].plot()
        for pipe in self.pipes[1:]:
            pipe.plot(ax=ax)
        return ax

    def volmdlr_primitives(self):
        pipes_primitives = []
        for pipe in self.pipes:
            if hasattr(pipe, 'volmdlr_primitives'):
                pipes_primitives.extend(pipe.volmdlr_primitives())
        return pipes_primitives

    @classmethod
    def dict_to_object(cls, dict_):
        points = [vm.Point3D.dict_to_object(d) for d in dict_['points']]
        pipes = [hyp.Pipe.dict_to_object(d) for d in dict_['pipes']]
        fluid = fluids.Fluid.dict_to_object(dict_['fluid'])
        boundary_conditions = []
        for d in dict_['boundary_conditions']:
            if d['type'] == 'pressure':
                PressureCondition
        circuit = cls(points, pipes, boundary_conditions, fluid)
        return circuit


class BoundaryCondition(dcc.DessiaObject):
    _standalone_in_db = False

    def __init__(self, point, value):
        self.points = [point]
        self.active_points = self.points
        self.value = value
        self.heat_exchange = False
        # self.system_matrix = self.system_matrix()
        self.n_equations = 1


class PressureCondition(BoundaryCondition):
    def __init__(self, point, value):
        BoundaryCondition.__init__(self, point, value)
        
    def __str__(self):
        return 'P={}'.format(round(self.value, 3))
        
    def system_matrix(self):
        system_matrix = npy.array([[1, 0]])
        return system_matrix

    def to_dict(self):
        d = copy(self.__dict__)
        d['points'] = [point.to_dict() for point in self.points]
        d['active_points'] = [point.to_dict() for point in self.active_points]
        d['system_matrix'] = self.system_matrix.tolist()
        d['type'] = 'pressure'
        return d

    @classmethod
    def dict_to_object(cls, dict_):
        point = [vm.Point3D.dict_to_object(d) for d in dict_['points']][0]
        value = dict_['value']
        condition = cls(point, value)
        return condition


class FlowCondition(dcc.DessiaObject):
    _standalone_in_db = True
    
    def __init__(self, point, value):
        BoundaryCondition.__init__(point, value)
        
    def __str__(self):
        return 'Q={}'.format(round(self.value, 3))

    def system_matrix(self):
        system_matrix = npy.array([[0, -1]])
        return system_matrix

    def to_dict(self):
        d = copy(self.__dict__)
        d['points'] = [point.to_dict() for point in self.points]
        d['active_points'] = [point.to_dict() for point in self.active_points]
        d['system_matrix'] = self.system_matrix.tolist()
        d['type'] = 'flow'
        return d

    @classmethod
    def dict_to_object(cls, dict_):
        point = [vm.Point3D.dict_to_object(d) for d in dict_['points']][0]
        value = dict_['value']
        condition = cls(point, value)
        return condition

class FluidicsResults(dcc.DessiaObject):
    
    def __init__(self, circuit, solution):
        self.circuit = circuit
        self.solution = solution
        system_data = circuit.settings()
        self.flows_dict = system_data[1]
        
    def to_thermal(self, fluid_input_points, fluid_output_points):
        """
        Converts hydraulic circuit to thermal
        """
        # Equivalence dictionaries initialising
        point2node = {}
        pipe2block = {}

        # Nodes/Blocks lists initialising
        nodes = []
        blocks = []
        wall_nodes = []

        # Counts initialising
        node_count = 0
        block_count = 0

        # Loop on every hydraulic pipe to create thermal equivalence
        for pipe in self.circuit.pipes:
            block_nodes = []
            for point in pipe.active_points:
                if point not in point2node.keys():
                    node = th.Node('h'+str(node_count))
                    node_count += 1
                    point2node[point] = node
                else:
                    node = point2node[point]
                block_nodes.append(node)
            if pipe.__class__ == hyp.JunctionPipe:
                # Junction pipe doesn't exchange heat
                flows_dict = [self.flows_dict[(pipe, point)] for point in pipe.active_points]
                flows = [self.solution[j] for j in flows_dict]
                input_indices = [index for index, flow in enumerate(flows) if flow >= 0]
                input_nodes = [node for index, node in enumerate(block_nodes) if index in input_indices]
                output_nodes = [node for index, node in enumerate(block_nodes) if index not in input_indices]
                input_flows = [flow for index, flow in enumerate(flows) if index in input_indices]
                block = th.Junction(input_nodes,
                                    output_nodes,
                                    input_flows,
                                    self.circuit.fluid,
                                    'b'+str(block_count))

            else:
                # Other pipes
                j = self.flows_dict[(pipe, pipe.points[0])]
                th_node = th.Node('t'+str(block_count))
                wall_nodes.append(th_node)
                nodes.append(th_node)
                block = th.ThermalPipe(block_nodes + [th_node],
                                       self.solution[j],
                                       self.circuit.fluid,
                                       'b'+str(block_count))

            block_count += 1
            pipe2block[pipe] = block
            for node in block.nodes:
                if node not in nodes:
                    nodes.append(node)
            blocks.append(block)

        interface_nodes = {'input': [point2node[fip] for fip in fluid_input_points],
                           'wall_nodes': wall_nodes,
                           'output': [point2node[fop] for fop in fluid_output_points]}

        thermal_circuit = th.Circuit(nodes, blocks)
        thermohydraulic_ciruit = th.ThermohydraulicCircuit(self.circuit,
                                                           thermal_circuit,
                                                           interface_nodes,
                                                           point2node,
                                                           pipe2block)
        return thermohydraulic_ciruit

    def plot(self, x3D=vm.X3D, y3D=vm.Y3D, position=False, color_map=None, max_width=0.025):
        """
        Displays the solution : [pressure,flow rate]
        Numbers are displayed with the number of digits
        position=True to display in circuit coordinates
        """
        dictionaries = self.circuit.settings()
        pressures_dict, flows_dict = dictionaries[:2]
        graph = self.circuit.graph

        points = [point for point in graph.nodes if graph.nodes[point]['node_type'] == 'point']
        pipes = [pipe for pipe in graph.nodes if graph.nodes[pipe]['node_type'] == 'pipe']
        pipe_points_couples = [(pipe, point) for pipe in pipes for point in pipe.active_points]
        flows = {}

        if position:
            # Points positions
            graph_pos = {point: (point.dot(x3D),point.dot(y3D))
                         for point in list(self.circuit.graph.nodes)\
                         if graph.nodes[point]['node_type'] == 'point'}

            # Pipes positions : centroid of points
            for pipe in pipes:
                # coords = [(point.x, point.y) for point in pipe.active_points]
                length = len(pipe.active_points)
                sum_x = sum([point.x for point in pipe.active_points])
                sum_y = sum([point.y for point in pipe.active_points])
                # sum_y = npy.sum(coords[:, 1])
                graph_pos[pipe] = (sum_x/length, sum_y/length)
        else:
            graph_pos = nx.kamada_kawai_layout(self.circuit.graph)

        fig = plt.figure()
        ax = fig.add_subplot(111)

        for ppc in pipe_points_couples:
            if ppc in flows_dict.keys():
                index = flows_dict[ppc]
                flows[ppc] = self.solution[index]
            else:
                edge_position = (graph_pos[ppc[0]], graph_pos[ppc[1]])
                x = [edge_position[0][0], edge_position[1][0]]
                y = [edge_position[0][1], edge_position[1][1]]
                ax.plot(x, y, 'k-')

        width = float(max_width/max(flows.values()))

        for ppc, flow in flows.items():
            add_link = False
            if flow < 0:
                # From block to node
                edge_position = (graph_pos[ppc[0]], graph_pos[ppc[1]])
            elif flow > 0:
                # From node to block
                edge_position = (graph_pos[ppc[1]], graph_pos[ppc[0]])
            else:
                add_link = True

            if add_link:
                edge_position = (graph_pos[ppc[0]], graph_pos[ppc[1]])
                x = [edge_position[0][0], edge_position[1][0]]
                y = [edge_position[0][1], edge_position[1][1]]
                ax.plot(x, y, 'k-')
            else:
                x = edge_position[0][0]
                y = edge_position[0][1]
                dx = edge_position[1][0] - x
                dy = edge_position[1][1] - y
                ax.add_patch(patches.FancyArrow(x, y, dx, dy, width=abs(width*flow),
                                                length_includes_head=True,
                                                head_width=abs(width*flow)*2,
                                                edgecolor='k',
                                                facecolor='gray'))

        press_values = [self.solution[pressures_dict[point]] for point in points]

        press_pts = nx.draw_networkx_nodes(graph,
                                           pos=graph_pos,
                                           nodelist=points,
                                           node_color=press_values,
                                           node_size=100,
                                           cmap=color_map)
        nx.draw_networkx_nodes(graph,
                               pos=graph_pos,
                               nodelist=pipes,
                               node_shape='s',
                               node_color='k',
                               node_size=50)

        ax.axis('equal')
        plt.colorbar(press_pts)

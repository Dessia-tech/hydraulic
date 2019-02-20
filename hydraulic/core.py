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
from copy import copy

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

class Circuit:
    """
    General class for 2D/3D circuits
    """
    def __init__(self, points, pipes, boundary_conditions, fluid):
        self.points = points
        self.pipes = pipes
        self.boundary_conditions = boundary_conditions
        self.fluid = fluid

        self.graph = self.GenerateGraph()
        # Creating thermal nodes out of circuit points

    def GenerateGraph(self):
        graph = nx.Graph()
        if self.points and self.pipes:
            active_points = {point for pipe in self.pipes\
                             for point in pipe.active_points}
            graph.add_nodes_from(active_points,
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

    def AddPoint(self, coordinates):
        point = vm.Point2D(coordinates)
        self.points.append(point)
        self.graph.add_node(point)

    def AddPipe(self, p1, p2, pipe_type, *infos):
        # Adds a pipe knowing the points linked, the type of pipe and the defining characterictics
        pipe = globals()[pipe_type](p1, p2, *infos)
        self.pipes.append(pipe)
        self.graph.add_edge(*pipe.points, object=pipe)

    def AddJunction(self, p1, ps, diameter, uniting=1):
        # uniting=1 if flows from ps to p1, =0 if from p1 to ps
        for p in ps:
            if uniting:
                pipe = hyp.JunctionPipe(p, p1, diameter)
            else:
                pipe = hyp.JunctionPipe(p1, p, diameter)
            self.pipes.append(pipe)
            self.graph.add_edge(*pipe.points, object=pipe)

    def Point2Pipes(self):
        point2pipes = {point : [] for point in self.points}
        for point in self.points:
            for pipe in self.pipes + self.boundary_conditions:
                if point in pipe.points:
                    point2pipes[point].append(pipe)
        return point2pipes

    def DrawGraph(self):
        graph_pos = nx.kamada_kawai_layout(self.graph)
        labels_dict = {node : str(node) for node in list(self.graph.nodes)}
        plt.figure("circuit")
        nx.draw_networkx(self.graph,
                         pos=graph_pos,
                         labels=labels_dict,
                         withlabels=True,
                         fontweight="bold")

    def ResolutionSettings(self):
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
           
        points2pipes = self.Point2Pipes()
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

    def SystemMatrix(self, constant):
        pressures_dict, flows_dict, equations_dict, neqs, nvars = self.ResolutionSettings()
        system_matrix = npy.zeros((neqs, nvars))
        points2pipes = self.Point2Pipes()

        # Write blocks equations
        for pipe in self.pipes + self.boundary_conditions:
            if pipe in self.boundary_conditions:
                block_system_matrix = pipe.SystemMatrix()
            else:
                block_system_matrix = pipe.SystemMatrix(constant)
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

    def SolveFluidics(self):
        """
        solves the hydraulic circuit for
        :param imposed: = "pressure" or "flow"
        pressure : imposed_values = deltaP = dictionary {point:pressure}
        flow : imposed_values = Qs = dictionary {point:flow rate}
               flow is set at input points
               for output points specify "out" as value
               output points will be set at pressure = low_p
        """
        dictionaries = self.ResolutionSettings()
        equations_dict = dictionaries[2]

        constant = 2/(self.fluid.rho*self.fluid.nu)
        system_matrix = self.SystemMatrix(constant)
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
            if Re > 2000: # Laminar condition
                turbulent = 1
                print("Warning, flow is not laminar in {} with Re = {}".format(pipe, Re))
        if not turbulent:
            print("Flow is laminar in the circuit")

    def Dict(self):
        d = {'points' : [point.Dict() for point in self.points]}
        d['pipes'] = [pipe.Dict() for pipe in self.pipes]
        d['boundary_conditions'] = [bc.Dict() for bc in self.boundary_conditions]
        d['fluid'] = self.fluid.Dict()
        return d

class Circuit2D(Circuit):
    def __init__(self, points, pipes, boundary_conditions, fluid):
        Circuit.__init__(self, points, pipes, boundary_conditions, fluid)

    def Draw(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        for pipe in self.pipes:
            pipe.Draw(ax)

    def Export1D(self, name="Generated_Circuit", path="", size=1., index=0):
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
            geo_txt += "Point({0}) = [{1[0]},{1[1]},0.,{2}];\n".format(i+1, point.position, size)
            pts_index[point] = i+1
        n_pipes = len(self.pipes)
        added_elems = 0
        physicals = [] # Index of lines for physical lines
        for j in range(n_pipes):
            pipe = self.pipes[j]
            [mess, phys] = pipe.Repr1D(j+1+added_elems, pts_index)
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
    3D circuit
    
    """
    
    def __init__(self, points, pipes, boundary_conditions, fluid):
        Circuit.__init__(self, points, pipes, boundary_conditions, fluid)

    def Draw(self, x3D=vm.x3D, y3D=vm.y3D, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        for pipe in self.pipes:
            pipe.Draw(x3D, y3D, ax)

    def CADModel(self):
        pipes_primitives = []
        for pipe in self.pipes:
            if hasattr(pipe, 'CADVolume'):
                pipes_primitives.append(pipe.CADVolume())
        model = vm.VolumeModel([('pipes', pipes_primitives)])
        return model

    def FreeCADExport(self, filename='An_unamed_circuit',
                      python_path='python',
                      path_lib_freecad='/usr/lib/freecad/lib',
                      export_types=('fcstd')):
        model = self.CADModel()
        return model.FreeCADExport(filename, python_path, path_lib_freecad, export_types)

    @classmethod
    def DictToObject(cls, dict_):
        points = [vm.Point3D.DictToObject(d) for d in dict_['points']]
        pipes = [hyp.Pipe.DictToObject(d) for d in dict_['pipes']]
        fluid = fluids.Fluid.DictToObject(dict_['fluid'])
        boundary_conditions = []
        for d in dict_['boundary_conditions']:
            if d['type'] == 'pressure':
                PressureCondition
        circuit = cls(points, pipes, boundary_conditions, fluid)
        return circuit


class BoundaryCondition:
    def __init__(self, point, value):
        self.points = [point]
        self.active_points = self.points
        self.value = value
        self.heat_exchange = False
        self.system_matrix = self.SystemMatrix()
        self.n_equations = 1

class PressureCondition(BoundaryCondition):
    def __init__(self, point, value):
        BoundaryCondition.__init__(self, point, value)
        
    def SystemMatrix(self):
        system_matrix = npy.array([[1, 0]])
        return system_matrix

    def Dict(self):
        d = copy(self.__dict__)
        d['points'] = [point.Dict() for point in self.points]
        d['active_points'] = [point.Dict() for point in self.active_points]
        d['system_matrix'] = self.system_matrix.tolist()
        d['type'] = 'pressure'
        return d

    @classmethod
    def DictToObject(cls, dict_):
        point = [vm.Point3D.DictToObject(d) for d in dict_['points']][0]
        value = dict_['value']
        condition = cls(point, value)
        return condition

class FlowCondition:
    def __init__(self, point, value):
        BoundaryCondition.__init__(point, value)

    def SystemMatrix(self):
        system_matrix = npy.array([[0, -1]])
        return system_matrix

    def Dict(self):
        d = copy(self.__dict__)
        d['points'] = [point.Dict() for point in self.points]
        d['active_points'] = [point.Dict() for point in self.active_points]
        d['system_matrix'] = self.system_matrix.tolist()
        d['type'] = 'flow'
        return d

    @classmethod
    def DictToObject(cls, dict_):
        point = [vm.Point3D.DictToObject(d) for d in dict_['points']][0]
        value = dict_['value']
        condition = cls(point, value)
        return condition

class FluidicsResults:
    def __init__(self, circuit, solution):
        self.circuit = circuit
        self.solution = solution
        system_data = circuit.ResolutionSettings()
        self.flows_dict = system_data[1]
        
    def ToThermal(self, fluid_input_points, fluid_output_points):
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
                input_nodes = [node for index, node in enumerate(block_nodes)\
                               if index in input_indices]
                output_nodes = [node for index, node in enumerate(block_nodes)\
                                if index not in input_indices]
                input_flows = [flow for index, flow in enumerate(flows)\
                               if index in input_indices]
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

        interface_nodes = {'input' : [point2node[fip] for fip in fluid_input_points],
                           'wall_nodes' : wall_nodes,
                           'output' : [point2node[fop] for fop in fluid_output_points]}

        thermal_circuit = th.Circuit(nodes, blocks)
        thermohydraulic_ciruit = th.ThermohydraulicCircuit(self.circuit,
                                                           thermal_circuit,
                                                           interface_nodes,
                                                           point2node,
                                                           pipe2block)
        return thermohydraulic_ciruit

    def DisplaySolution(self, x3D=vm.x3D, y3D=vm.y3D,
                        position=False,
                        color_map=None,
                        max_width=0.025):
        """
        Displays the solution : [pressure,flow rate]
        Numbers are displayed with the number of digits
        position=True to display in circuit coordinates
        """
        dictionaries = self.circuit.ResolutionSettings()
        pressures_dict, flows_dict = dictionaries[:2]
        graph = self.circuit.graph

        points = [point for point in graph.nodes if graph.nodes[point]['node_type'] == 'point']
        pipes = [pipe for pipe in graph.nodes if graph.nodes[pipe]['node_type'] == 'pipe']
        pipe_points_couples = [(pipe, point) for pipe in pipes\
                               for point in pipe.active_points]
        flows = {}

        if position:
            # Points positions
            graph_pos = {point : npy.array((npy.dot(point.vector, x3D.vector),
                                            npy.dot(point.vector, y3D.vector)))\
                         for point in list(self.circuit.graph.nodes)\
                         if graph.nodes[point]['node_type'] == 'point'}

            # Pipes positions : centroid of points
            for pipe in pipes:
                coords = npy.asarray([point.vector for point in pipe.active_points])
                length = coords.shape[0]
                sum_x = npy.sum(coords[:, 0])
                sum_y = npy.sum(coords[:, 1])
                graph_pos[pipe] = npy.array((sum_x/length, sum_y/length))
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
    
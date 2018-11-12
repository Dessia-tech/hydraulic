#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Hydraulic circuit definition
"""

import math
import numpy as npy
import networkx as nx
import matplotlib.pyplot as plt
from scipy.linalg import solve
import hydraulic.pipes as hy_pipes
import hydraulic.thermal as thermal
import volmdlr as vm

# Definition of equivalent L/D values
# lists in increasing order
ben_angle = [0, 45, 90, 180, 360]
ben_LD = [0, 16, 30, 60, 120]
# bend radius ?
mbd_angle = [0, 45, 90]
mbd_LD = [0, 15, 58]
enl_rap = [1/4, 1/2, 3/4, 1]
enl_LD = [30, 20, 6.5, 0]
ctr_rap = [1, 4/3, 2, 4]
ctr_LD = [0, 6.5, 11, 15]

class ThermalNode:
    def __init__(self, name=''):
        self.name = name

class ThermalResistance:
    def __init__(self, resistivity):
        self.resistivity = resistivity

class Circuit:
    """
    general class for 2D circuits
    """
    def __init__(self, points, pipes, fluid):
        self.points = points
        self.pipes = pipes
        self.fluid = fluid

        self.graph = self.GenerateGraph()
        # Creating thermal nodes out of circuit points

    def GenerateGraph(self):
        graph = nx.Graph()
        if self.points:
            graph.add_nodes_from(self.points)
            for pipe in self.pipes:
                graph.add_edge(*pipe.points, object=pipe)
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
                pipe = hy_pipes.JunctionPipe(p, p1, diameter)
            else:
                pipe = hy_pipes.JunctionPipe(p1, p, diameter)
            self.pipes.append(pipe)
            self.graph.add_edge(*pipe.points, object=pipe)

    def DrawGraph(self):
        graph_pos = nx.kamada_kawai_layout(self.graph)
        labels_dict = {node : str(node) for node in list(self.graph.nodes)}
        dic_edges = nx.get_edge_attributes(self.graph, 'object')
        labels_edge_dict = {edge : dic_edges[edge] for edge in list(self.graph.edges)}
        plt.figure("circuit")
        nx.draw_networkx(self.graph,
                         pos=graph_pos,
                         labels=labels_dict,
                         withlabels=True,
                         fontweight="bold")
        nx.draw_networkx_edge_labels(self.graph,
                                     pos=graph_pos,
                                     edge_labels=labels_edge_dict,
                                     withlabels=True,
                                     font_size=7)

    def SolveFluidics(self, imposed, imposed_values, low_p=0):
        """
        solves the hydraulic circuit for
        :param imposed: = "pressure" or "flow"
        pressure : imposed_values = deltaP = dictionary {point:pressure}
        flow : imposed_values = Qs = dictionary {point:flow rate}
               flow is set at input points
               for output points specify "out" as value
               output points will be set at pressure = low_p
        """
        if imposed == "pressure":
            delta_P = imposed_values
            Qs = dict()
            l_index_q = []
        elif imposed == "flow":
            Qs = imposed_values
            delta_P = dict()
            l_index_q_temp = [self.points.index(point) for point in list(Qs.keys())]
            l_index_q = []
            for i in l_index_q_temp:
                if  Qs[self.points[i]] == "out":
                    delta_P[self.points[i]] = low_p
                else:
                    l_index_q.append(i)
        else:
            print("Warning, must specify if flow or pressure is imposed")
            return ()
        n = len(self.points)
        m = len(self.pipes)
        M_behavior = npy.zeros([m, n])
        S_nodes = npy.zeros([n, m])
        eps = 1e-3 # Small value for approximation of deltaP=0
        # A too small eps value may result in Ill-conditioned matrix

        for j in range(m):
            pipe = self.pipes[j]
            i1 = self.points.index(pipe.points[0])
            i2 = self.points.index(pipe.points[1])
            if pipe.fQ:
                M_behavior[j, i1] = 2/(self.fluid.rho*self.fluid.nu*pipe.fQ)
                M_behavior[j, i2] = -2/(self.fluid.rho*self.fluid.nu*pipe.fQ)
            else:
                M_behavior[j, i1] = 2/(self.fluid.rho*self.fluid.nu*eps)
                M_behavior[j, i2] = -2/(self.fluid.rho*self.fluid.nu*eps)
            S_nodes[i1, j] = -1
            S_nodes[i2, j] = 1

        A_circuit = npy.dot(S_nodes, M_behavior)
        M_circuit = A_circuit.copy()
        l_index_p = [self.points.index(point) for point in list(delta_P.keys())]
        l_index_p.sort()
        n_del = 0
        for i in l_index_p:
            M_circuit = npy.delete(M_circuit, i-n_del, axis=0)
            n_del += 1
        vect_delta_P = npy.zeros([n-n_del,])
        vect_imp_Q = npy.zeros([n-n_del,])
        n_del_col = 0
        l_index = l_index_p+l_index_q
        l_index.sort()
        for i in l_index:
            if i in l_index_p:
                vect_delta_P += M_circuit[:, i-n_del_col]*(-1)*delta_P[self.points[i]]
                M_circuit = npy.delete(M_circuit, i-n_del_col, axis=1)
                n_del_col += 1
            else:
                vect_imp_Q[i-n_del_col] = -Qs[self.points[i]]
        vect_imp = vect_delta_P + vect_imp_Q
        vect_P = list(solve(M_circuit, vect_imp))
        for i in l_index_p:
            vect_P.insert(i, delta_P[self.points[i]])

        vect_Q = npy.dot(M_behavior, vect_P)
        # Add solution to each element
#        for i in range(0,n):
#            pressures[self.points[i]] = vect_P[i]
##            self.points[i].press=vect_P[i]
#        for i in range(0,m):
#            flows[self.points[i]] = vect_Q[i]
##            self.pipes[i].flow=vect_Q[i]
        result = FluidicsResults(vect_P, vect_Q, self, imposed, imposed_values, low_p)
        return result

    def VerifyQ(self):
        """
        Check if the flow rate is laminar inside the circuit
        """
        turbulent = 0
        for i, pipe in enumerate(self.pipes):
            Re = 2*pipe.flow/(self.fluid.nu*math.pi*pipe.radius)
            if Re > 2000: # Laminar condition
                turbulent = 1
                print("Warning, flow is not laminar in {} with Re = {}".format(pipe, Re))
        if not turbulent:
            print("Flow is laminar in the circuit")

    def ToThermal(self, resistors):
        cooling_plates = []
        for i, pipe in enumerate(self.pipes):
            cooling_plates.append(thermal.CoolingPlatePartFromPipe(pipe, resistors[i], self.fluid))
        return cooling_plates

class Circuit2D(Circuit):
    def __init__(self, points, pipes, fluid):
        Circuit.__init__(self, points, pipes, fluid)

    def Draw(self):
        fig, ax = plt.subplots()
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
            return(pts_index, n_pipes+added_elems, physicals)

    def FindLimits(self):
        # Returns the nodes which are defined as the end of a pipe
        points_lim = []
        for point in self.points:
            if len(self.graph.adj[point]) == 1:
                points_lim.append(point)
        return points_lim

class Circuit3D(Circuit):
    def __init__(self, points, pipes, fluid):
        Circuit.__init__(self, points, pipes, fluid)

    def Draw(self, x3D, y3D):
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        for pipe in self.pipes:
            pipe.Draw(x3D, y3D, ax)

    def CADModel(self):
        pipes_primitives = []
        for pipe in self.pipes:
            pipes_primitives.append(pipe.CADVolume())
        model = vm.VolumeModel([('pipes', pipes_primitives)])
        return model

    def FreeCADExport(self, filename='An_unamed_circuit',
                      python_path='python',
                      path_lib_freecad='/usr/lib/freecad/lib',
                      export_types=('fcstd')):
        model = self.CADModel()
        return model.FreeCADExport(filename, python_path, path_lib_freecad, export_types)

class FluidicsResults:
    def __init__(self, vect_p, vect_q,
                 circuit,
                 imposed, imposed_values, low_p):
        self.vect_p = vect_p
        self.vect_q = vect_q
        self.circuit = circuit
        self.imposed = imposed
        self.imposed_values = imposed_values
        self.low_p = low_p
        self.dict_press = {point : vect_p[i] for i, point in enumerate(circuit.points)}
        self.dict_flows = {pipe : vect_q[i] for i, pipe in enumerate(circuit.pipes)}

    def DisplaySolution(self, digits=2,
                        position=False,
                        color_map=None,
                        values=False,
                        max_width=7):
        """
        Displays the solution : [pressure,flow rate]
        Numbers are displayed with the number of digits
        position=True to display in circuit coordinates
        """
        if position:
            graph_pos = {node:node.vector for node in list(self.circuit.graph.nodes)}
            plt.figure("Solution in circuit coordinates")
        else:
            graph_pos = nx.kamada_kawai_layout(self.circuit.graph)
            plt.figure("Solution")

        if values:
            # Format Data
            dict_press = {point : "%.{}E".format(digits) % (pressure)\
                          for point, pressure in self.dict_press}
            dict_flows = {pipe : "%.{}E".format(digits) % (flow)\
                          for pipe, flow in self.dict_flows}

            labels_dict = {node : dict_press[node]\
                           for node in list(self.circuit.graph.nodes)}
            dict_edges = nx.get_edge_attributes(self.circuit.graph, 'object')
            labels_edge_dict = {edge : dict_flows[dict_edges[edge]]\
                                for edge in list(self.circuit.graph.edges)}

            # Display
            nx.draw_networkx(self.circuit.graph,
                             pos=graph_pos,
                             labels=labels_dict,
                             withlabels=True,
                             node_color='c')
            nx.draw_networkx_edge_labels(self.circuit.graph,
                                         pos=graph_pos,
                                         edge_labels=labels_edge_dict,
                                         withlabels=True,
                                         font_size=7)
        else:
            width = float(max_width/max(self.dict_flows.values()))
            press_color = [self.dict_press[node]\
                           for node in list(self.circuit.graph.nodes)]
            dic_edges = nx.get_edge_attributes(self.circuit.graph, 'object')
            flow_width = [width*self.dict_flows[dic_edges[edge]]\
                          for edge in list(self.circuit.graph.edges)]

            # Display
            press_pts = nx.draw_networkx_nodes(self.circuit.graph,
                                               pos=graph_pos,
                                               node_color=press_color,
                                               cmap=color_map)
            nx.draw_networkx_edges(self.circuit.graph,
                                   pos=graph_pos,
                                   width=0.5,
                                   style='dashed')
            nx.draw_networkx_edges(self.circuit.graph,
                                   pos=graph_pos,
                                   width=flow_width,
                                   edge_color='k')
            plt.colorbar(press_pts)
    
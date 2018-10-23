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

# Definition of equivalent L/D values
# lists in increasing order
ben_angle = [0,45,90,180,360]
ben_LD = [0,16,30,60,120]
# bend radius ?
mbd_angle = [0,45,90]
mbd_LD = [0,15,58]
enl_rap = [1/4,1/2,3/4,1]
enl_LD = [30,20,6.5,0]
ctr_rap = [1,4/3,2,4]
ctr_LD = [0,6.5,11,15]

from volmdlr import Point2D, Arc2D, LineSegment2D, Point3D, Arc3D, LineSegment3D

class Circuit:
    """
    general class for 2D circuits
    """
    def __init__(self,points ,pipes, fluid):
        self.points = points
        self.pipes = pipes
        self.fluid = fluid
        
        self.graph = self.GenerateGraph()
        
    def GenerateGraph(self):
        graph=nx.Graph()
        if self.points :
            graph.add_nodes_from(self.points)
            for pipe in self.pipes:
                graph.add_edge(*pipe.points, object=pipe)
        return(graph)
        
    def AddPoint(self,coordinates):
        point = Point2D(coordinates)
        self.points.append(point)
        self.graph.add_node(point)
        
    def AddPipe(self,p1,p2,pipe_type,*infos):
        # adds a pipe knowing the points linked, the type of pipe and the defining characterictics
        pipe=globals()[pipe_type](p1,p2,*infos)
        #print(pipe)
        self.pipes.append(pipe)
        self.graph.add_edge(*pipe.points,object=pipe)
        
    def AddJunction(self,p1,ps,diameter,uniting=1):
        # uniting=1 if flows from ps to p1, =0 if from p1 to ps
        for p in ps:
            if uniting:
                pipe=JunctionPipe(p,p1,diameter)
            else:
                pipe=JunctionPipe(p1,p,diameter)
            self.pipes.append(pipe)
            self.graph.add_edge(*pipe.points,object=pipe)
            
    def DrawGraph(self):
        graph_pos=nx.kamada_kawai_layout(self.graph)
        labels_dict={node:str(node) for node in list(self.graph.nodes)}
        dic_edges=nx.get_edge_attributes(self.graph,'object')
        labels_edge_dict={edge:dic_edges[edge] for edge in list(self.graph.edges)}
        plt.figure("circuit")
        nx.draw_networkx(self.graph,pos=graph_pos,labels=labels_dict,withlabels=True,fontweight="bold")
        nx.draw_networkx_edge_labels(self.graph,pos=graph_pos,edge_labels=labels_edge_dict,withlabels=True,font_size=7)
    

        
    def Solve(self,imposed, imposed_values,low_p=0):
        """
        solves the hydraulic circuit for 
        :param imposed: = "pressure" or "flow"
        pressure : imposed_values = deltaP = dictionary {point:pressure}
        flow : imposed_values = Qs = dictionary {point:flow rate}
               flow is set at input points
               for output points specify "out" as value
               output points will be set at pressure = low_p
        """
        if imposed=="pressure":
            delta_P=imposed_values
            Qs=dict()
            l_index_q=[]
        elif imposed=="flow":
            Qs=imposed_values
            delta_P=dict()
            l_index_q_temp=[self.points.index(pt) for pt in list(Qs.keys())]
            l_index_q=[]
            for i in l_index_q_temp:
                if  Qs[self.points[i]]=="out":
                    delta_P[self.points[i]]=low_p
                else:
                    l_index_q.append(i)
        else:
            print("Warning, must specify if flow or pressure is imposed")
            return()
        n=len(self.points)
        m=len(self.pipes)
        M_behavior=npy.zeros([m,n])
        S_nodes=npy.zeros([n,m])
        eps=1e-3 #small value for approximation of deltaP=0
        # eps too small may result in Ill-conditioned matrix
        for j in range(m):
            pipe=self.pipes[j]
            i1=self.points.index(pipe.points[0])
            i2=self.points.index(pipe.points[1])
            if pipe.fQ:
                M_behavior[j,i1]=2/(self.fluid.rho*self.fluid.nu*pipe.fQ)
                M_behavior[j,i2]=-2/(self.fluid.rho*self.fluid.nu*pipe.fQ)
            else:
                M_behavior[j,i1]=2/(self.fluid.rho*self.fluid.nu*eps)
                M_behavior[j,i2]=-2/(self.fluid.rho*self.fluid.nu*eps)
            S_nodes[i1,j]=-1
            S_nodes[i2,j]=1

        A_circuit=npy.dot(S_nodes,M_behavior)
        M_circuit=A_circuit.copy()
        l_index_p=[self.points.index(pt) for pt in list(delta_P.keys())]
        l_index_p.sort()
        n_del=0
        for i in l_index_p:
            M_circuit=npy.delete(M_circuit,i-n_del,axis=0)
            n_del+=1
        vect_delta_P=npy.zeros([n-n_del,])
        vect_imp_Q=npy.zeros([n-n_del,])
        n_del_col=0
        l_index=l_index_p+l_index_q
        l_index.sort()
        for i in l_index:
            if i in l_index_p:
                vect_delta_P+=M_circuit[:,i-n_del_col]*-1*delta_P[self.points[i]]
                M_circuit=npy.delete(M_circuit,i-n_del_col,axis=1)
                n_del_col+=1
            else:
                vect_imp_Q[i-n_del_col]=-Qs[self.points[i]]
        vect_imp=vect_delta_P+vect_imp_Q
        vect_P=list(solve(M_circuit,vect_imp))
        for i in l_index_p:
            vect_P.insert(i,delta_P[self.points[i]])
        #vect_P=npy.reshape(vect_P,[n-2,1])
        vect_Q=npy.dot(M_behavior,vect_P)
        # add solution to each element
        for i in range(0,n):
            self.points[i].press=vect_P[i]
        for i in range(0,m):
            self.pipes[i].flow=vect_Q[i]
        return(vect_P,vect_Q)
        
        
    def VerifyQ(self):
        """
        Check if the flow rate is laminar inside the circuit
        """
        m=len(self.pipes)
        turbulent=0
        for i in range(m):
            Re=2*self.pipes[i].flow/(self.nu*math.pi*self.pipes[i].radius)
            if Re>2000: # laminar condition
                turbulent=1
                print("Warning, flow is not laminar in {} with Re = {}".format(self.pipes[i],Re))
        if not turbulent:
            print("Flow is laminar in the circuit")
            
    def DisplaySolution(self,digits=2, position=0, color_map=None, values=0, max_width=7):
        """
        displays the solution=[pressure,flow rate]
        numbers are displayed with the numeber of digits
        position=1 for display in circuit coordinates
        """
        
        dic_press = dict()
        dic_flow = dict()
        if position:
            graph_pos={node:node.vector for node in list(self.graph.nodes)}
            plt.figure("solution in coordinates")
        else:
            graph_pos=nx.kamada_kawai_layout(self.graph)
            plt.figure("solution")
        if values:
            for pt in self.points:
                dic_press[pt]="%.{}E".format(digits) % (pt.press)
            for pipe in self.pipes:
                dic_flow[pipe]="%.{}E".format(digits) % (pipe.flow)
            labels_dict={node:dic_press[node] for node in list(self.graph.nodes)}
            dic_edges=nx.get_edge_attributes(self.graph,'object')
            labels_edge_dict={edge:dic_flow[dic_edges[edge]] for edge in list(self.graph.edges)}
            nx.draw_networkx(self.graph,pos=graph_pos,labels=labels_dict,withlabels=True,node_color='c')
            nx.draw_networkx_edge_labels(self.graph,pos=graph_pos,edge_labels=labels_edge_dict,withlabels=True,font_size=7)
        else:
            for pt in self.points:
                dic_press[pt]=pt.press
            for pipe in self.pipes:
                dic_flow[pipe]=pipe.flow
            width=float(max_width/max(dic_flow.values()))
            press_color=[dic_press[node] for node in list(self.graph.nodes)]
            dic_edges=nx.get_edge_attributes(self.graph,'object')
            flow_width=[width*dic_flow[dic_edges[edge]] for edge in list(self.graph.edges)]
            press_pts=nx.draw_networkx_nodes(self.graph,pos=graph_pos,node_color=press_color,cmap=color_map)
            nx.draw_networkx_edges(self.graph,pos=graph_pos,width=0.5,style='dashed')
            nx.draw_networkx_edges(self.graph,pos=graph_pos,width=flow_width,edge_color='k')
            plt.colorbar(press_pts)
            

    
class Circuit2D(Circuit):
    def __init__(self, points, pipes, fluid):
        Circuit.__init__(self, points, pipes, fluid)
        
    def Draw(self):
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        for pipe in self.pipes:
            pipe.Draw(ax)
        
    def Export1D(self,name="Generated_Circuit",path="",size=1.,index=0):
        """
        exports the circuit as a 1D geo file for gmsh
        size defines the size of the mesh
        """
        file_name=path+name+".geo"
        file=open(file_name,"w")
        n_pts=len(self.points)
        pts_index=dict()
        geo_txt=""
        for i in range(n_pts):
            pt=self.points[i]
            geo_txt+="Point({0}) = [{1[0]},{1[1]},0.,{2}];\n".format(i+1,pt.position,size)
            pts_index[pt]=i+1
        n_pipes=len(self.pipes)
        added_elems=0
        physicals=[] # index of lines for physical lines 
        for j in range(n_pipes):
            pipe=self.pipes[j]
            [mess,phys]=pipe.Repr1D(j+1+added_elems,pts_index)
            physicals.append(phys)
            geo_txt+=mess
            added_elems+=len(phys)-1
        geo_txt=geo_txt.replace("[","{")
        geo_txt=geo_txt.replace("]","}")
        file.write(geo_txt)
        file.close()
        if index:
            return(pts_index,n_pipes+added_elems,physicals)
    def FindLimits(self):
        #returns the nodes which are defined as the end of a pipe
        pts_lim=[]
        for pt in self.points:
            if len(self.graph.adj[pt])==1:
                pts_lim.append(pt)
        return pts_lim
    
    
class Circuit3D(Circuit):
    def __init__(self, points, pipes, fluid):
        Circuit.__init__(self, points, pipes, fluid)
        
    def Draw(self, x3D, y3D):
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        for pipe in self.pipes:
            print(pipe)
            pipe.Draw(x3D, y3D, ax)



class StraightPipe:
    """
    Abstract Class
    Straight pipes linking 2 points
    """
    def __init__(self, p1, p2, d):
        self.points = [p1, p2]
        self.radius = d/2
        self.surf = math.pi*self.radius**2
        self.length = p1.PointDistance(p2)
        self.fQ = 16*self.length/(math.pi*self.radius**4)
        
    def __repr__(self):
        return("{}-str-{}".format(self.points[0],self.points[1]))
        
    def __str__(self):
        return("{}-str-{}".format(self.points[0],self.points[1]))
        
    def Repr1D(self,j,points_index):
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])


class StraightPipe2D:
    """
    Straight pipes linking 2 2D points
    """
    def __init__(self, p1, p2, d):
        StraightPipe.__init__(self, p1, p2, d)

    def Draw(self, ax):
        l = LineSegment2D(*self.points)
        l.MPLPlot(ax)

        
class StraightPipe3D:
    """
    Straight pipes linking 2 3D points
    """
    def __init__(self, p1, p2, d):
        StraightPipe.__init__(self, p1, p2, d)

    def Draw(self, x3D, y3D, ax):
        l = LineSegment3D(*self.points)
        l.MPLPlot2D(x3D, y3D, ax)    
        
class SingularPipe:
    """
    other type of pipes linking 2 points
    """
    def __init__(self,p1, p2, form):
        self.points = [p1, p2]
        self.type=form
        
    def __repr__(self):
        return("{}-{}-{}".format(self.points[0],self.type,self.points[1]))
        
    def __str__(self):
        return("{}-{}-{}".format(self.points[0],self.type,self.points[1]))
    
class Bend(SingularPipe):
    """
    Bend between two points, defined by third, interior point
    """
    def __init__(self, start_point, interior_point, end_point, arc, diameter):
        SingularPipe.__init__(self, start_point, end_point, 'ben')
        self.start_point = start_point
        self.interior_point = interior_point
        self.end_point = end_point
        self.radius = 0.5 * diameter
#        [self.turn_radius,self.turn_angle,self.center] = ArcParam(p1.vector, p2.vector, coord_p3)
        self.arc = arc
        self.turn_radius = self.arc.radius
        self.turn_angle = self.arc.angle
        self.section = math.pi*self.radius**2
        self.length = self.turn_radius*self.turn_angle
        length_d = GetEquivalent(abs(self.turn_angle*180/math.pi) ,ben_LD, ben_angle)
        self.fQ=16*2*length_d/(math.pi*self.radius**3)
                
        

        
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        n=max(points_index.values())+1
        circle="Point({0}) = [{1[0]},{1[1]},0.,1.];\n".format(n+1,self.center)
        points_index["center of {}".format(str(self))]=n+1
        circle+="Point({0}) = [{1[0]},{1[1]},0.,1.];\n".format(n+2,self.third_point)
        points_index["intermediate point of {}".format(str(self))]=n+2
        circle+="Circle({}) = [{},{},{}];\n".format(j,points_index[self.points[0]],n+1,n+2)
        circle+="Circle({}) = [{},{},{}];\n".format(j+1,n+2,n+1,points_index[self.points[1]])
        phys=[j,j+1]
#        if self.turn_angle<math.pi:
#            circle+="Circle({}) = [{},{},{}];\n".format(j,points_index[self.points[0]],n+1,points_index[self.points[1]])
#            phys=[j]
#        else:
#            circle+="Point({0}) = [{1[0]},{1[1]},0.,1.];\n".format(n+2,self.third_point)
#            points_index["intermediate point of {}".format(str(self))]=n+2
#            circle+="Circle({}) = [{},{},{}];\n".format(j,points_index[self.points[0]],n+1,n+2)
#            circle+="Circle({}) = [{},{},{}];\n".format(j+1,n+2,n+1,points_index[self.points[1]])
#            phys=[j,j+1]
        return(circle,phys)
        
class Bend2D(Bend):
    """
    Bend between two points, defined by third, interior point
    """
    def __init__(self, start_point, interior_point, end_point, diameter):
        arc = Arc2D(start_point, interior_point, end_point)
        Bend.__init__(self, start_point, interior_point, end_point, arc, diameter)
        
    def Draw(self, ax):
        self.arc.MPLPlot(ax)

class Bend3D(Bend):
    """
    Bend between two points, defined by third, interior point
    """
    def __init__(self, start_point, interior_point, end_point, diameter):
        arc = Arc3D(start_point, interior_point, end_point)
        Bend.__init__(self, start_point, interior_point, end_point, arc, diameter)

    def Draw(self, x3D, y3D, ax):
        self.arc.MPLPlot2D(x3D, y3D, ax)    


class MitterBend(SingularPipe):
    """
    mitter bend 
    (for one point, specify 2 points with same coordinates)
    """
    def __init__(self,p1, p2,diameter,angle):
        SingularPipe.__init__(self,p1, p2,'mbd')
        self.radius=diameter/2
        self.surf=math.pi*self.radius**2
        self.length = p1.PointDistance(p2)
        length_d = GetEquivalent(angle*180/math.pi,mbd_LD,mbd_angle)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)
        
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])
        
class Enlargement(SingularPipe):
    """
    Enlargement from diameter 1 at point 1 to diameter 2 at point 2
    for gradual enlargement, specify the angle
    """
    def __init__(self,p1,p2,d1,d2,angle=None):
        SingularPipe.__init__(self, p1, p2, 'enl')
        if angle :
            s_angle = math.sin(angle)
        else :
            s_angle = 0
        self.radius=d1/2
        self.surf=math.pi*(d1/2)**2
        self.length = p1.PointDistance(p2)
        length_d=GetEquivalent(d1/d2,enl_LD,enl_rap)
        self.fQ=16*2*length_d/(math.pi*self.radius**3)
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])
        
class Contraction(SingularPipe):
    """
    Contraction from diameter 1 at point 1 to diameter 2 at point 2
    for gradual ontraction, specify the angle
    """
    def __init__(self,p1,p2,d1,d2,angle=None):
        SingularPipe.__init__(self,p1,p2,'ctr')
     
        if angle :
            s_angle=math.sin(angle)
        else :
            s_angle=0
        self.radius=d1/2
        self.surf=math.pi*(d1/2)**2
        self.length= p1.PointDistance(p2)
        length_d = GetEquivalent(d1/d2,ctr_LD,ctr_rap)
        self.fQ = 16*2*length_d/(math.pi*self.radius**3)
    
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])
        
class UserDefined(SingularPipe):
    """
    Singular Pipe defined by user with a known parameter and different 
    L/D equivalent ratios for different values of the parameter.
    """
    def __init__(self, p1, p2, diameter, param_values, LD_values, *params):
        SingularPipe.__init__(self, p1, p2, 'usr')
        self.radius=diameter/2
        self.surf=math.pi*self.radius**2
        self.length = p1.PointDistance(p2)
        length_d=GetEquivalent(param_values,LD_values,*params)
        self.fQ=16*2*length_d/(math.pi*self.radius**3)
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])
        
class JunctionPipe:
    """
    add pressure drop values
    junction linking 1 pipe to 2+ pipes
    """
    def __init__(self, p1, p2, diameter, value_LD=0):
        self.points = [p1, p2]
        self.radius = diameter/2
        self.surf = math.pi*self.radius**2
        self.length = p1.PointDistance(p2)
        length_d=self.length/diameter # to be replaced by values for junction pipes
        self.fQ=16*2*length_d/(math.pi*self.radius**3)
        
    def __repr__(self):
        return("{}-jun-{}".format(self.points[0],self.points[1]))
        
    def __str__(self):
        return("{}-jun-{}".format(self.points[0],self.points[1]))
        
    def Draw(self, ax):
        LineSegment2D(*self.points).MPLPlot(ax)
        
    def Repr1D(self,j,points_index):
        # returns the 1D representation of the pipe for gmsh
        line="Line({}) = [{},{}];\n".format(j,points_index[self.points[0]],points_index[self.points[1]])
        return(line,[j])
 
vm_equivalences = {LineSegment3D: (StraightPipe3D, (('points', 0), ('points', 1))),
                   Arc3D: (Bend3D, (('start', None), ('interior', None), ('end', None)))}

def PipesFromVolmdlrPrimitives(primitive, d):
    hy_class, args = vm_equivalences[primitive.__class__]
    args2 = []
    for a, index in args:
        if index is None:
            args2.append(getattr(primitive, a))
        else:
            args2.append(getattr(primitive, a)[index])
    print(args)
    return hy_class(*args2, d)
    
def GetEquivalent(pt_coord, known_val, *grid_val):
    # value for pt_coord using linear interpolation from discrete values inside a grid
    # y_val can be a list or numpy array
    # in case known_val is numpy array dim 2 then grid_val is 2 lists for lines and columns
    dim=len(npy.shape(known_val))
    if dim==1:
        Lx=grid_val[0]
        x=pt_coord
        y_val=known_val
        n,i=len(Lx),0
        while i<n and Lx[i]<x:
            i+=1
        if i==n:
            y=y_val[-1]
        else:
            y=y_val[i-1]+(x-Lx[i-1])*(y_val[i]-y_val[i-1])/(Lx[i]-Lx[i-1])
        return(y)
    elif dim==2:
        [x,y]=pt_coord
        [Lx,Ly]=grid_val
        nx,i=len(Lx),0
        while i<nx and Lx[i]<x:
            i+=1
        ny,j=len(Ly),0
        while j<ny and Ly[i]<y:
            j+=1
        i,j=i-1,j-1
        if i==nx-1 and j==ny-1:
            z=known_val[i,j]
        elif i==nx-1:
            z=known_val[i,j]+(y-Ly[j])*(known_val[i,j+1]-known_val[i,j])/(Ly[j+1]-Ly[j])
        elif j==ny-1:
            z=known_val[i,j]+(x-Lx[i])*(known_val[i+1,j]-known_val[i,j])/(Lx[i+1]-Lx[i])
        else:
            z=known_val[i,j]
            # interpolation on x
            z+=(known_val[i+1,j]-known_val[i,j])*(x-Lx[i])/(Lx[i+1]-Lx[i])
            # interpolation on y
            z+=(known_val[i,j+1]-known_val[i,j])*(y-Ly[j])/(Ly[j+1]-Ly[i])
            # interpolation on xy
            z+=(known_val[i+1,j+1]+known_val[i,j]-known_val[i,j+1]-known_val[i+1,j])*((y-Ly[j])/(Ly[j+1]-Ly[i]))*((x-Lx[i])/(Lx[i+1]-Lx[i]))
        return(z)
    else:
        print("Function GetEquivalent doesn't support dimension {}".format(dim))
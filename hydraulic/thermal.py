#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 17:54:15 2018

@author: jezequel
"""

import hydraulic as hy
import hydraulic.pipes as pipes
from hydraulic.fluids import water
import volmdlr as vm
import matplotlib.pyplot as plt

class ThermalResistor:
    def __init__(self, thickness, thermal_conductivity, area):
        self.thickness = thickness
        self.thermal_conductivity = thermal_conductivity
        self.area = area
        self.resistance = thickness/(thermal_conductivity*area)

    def Flow(self, delta_temperature):
        flow = delta_temperature/self.resistance
        return flow

    def DeltaTemperature(self, flow):
        delta_temperature = self.resistance*flow
        return delta_temperature

class CoolingPlatePart:
    def __init__(self, thermal_resistor, pipe, fluid):
        self.thermal_resistor = thermal_resistor
        self.pipe = pipe
        self.fluid = fluid

    def SolveThermic(self, T1, Tc, heat_flow):
        fluid_flow = self.pipe.flow
        rho = self.fluid.rho
        c = self.fluid.heat_capacity
        rth = self.thermal_resistor.resistance
        
        if heat_flow is not None and Tc is None: # Heat flow imposed
            delta_temperature = self.thermal_resistor.DeltaTemperature(heat_flow)
            T2 = delta_temperature*rth/(rho*c*fluid_flow) + T1
            Tp = (T2 + T1)/2.
            Tc = delta_temperature + Tp
        elif heat_flow is None and Tc is not None: # Cell temperature imposed
            T2 = (rth/(rho*c*fluid_flow + rth/2.))*Tc + T1
            Tp = (T2 + T1)/2
            delta_temperature = Tc - Tp
            heat_flow = self.thermal_resistor.Flow(delta_temperature)
        elif heat_flow is not None and Tc is not None:
            print("Warning ! Cannot impose both temperature and heat_flow")
            return False
        else:
            print("Warning ! Must impose temperature or heat_flow")
            return False
        return T2, heat_flow
    
class CoolingPlateJunction:
    def __init__(self, junction, fluid):
        self.junction = junction
        self.fluid = fluid
        
    def SolveThermic(self, in_temperatures, in_flows, out_flows):
        out_temperature = sum([t*f for t, f in zip(in_temperatures, in_flows)])/sum(in_flows)
        return out_temperature

class CoolingPlate:
    def __init__(self, circuit, cooling_plate_parts):
        self.circuit = circuit
        self.cooling_plate_parts = cooling_plate_parts

    def SolveFluidics(self, imposed, imposed_values):
        self.circuit.SolveFluidics(imposed, imposed_values)
    
    def SolveThermic(self, dp):
        for icp, cooling_plate_part in enumerate(self.cooling_plate_parts):
            a = cooling_plate_part.SolveThermic()
        return a
    
def CoolingPlatePartFromPipe(pipe, thermal_resistor, fluid):
    if issubclass(pipe.__class__, hy.pipes.StraightPipe)\
    or issubclass(pipe.__class__, hy.pipes.Bend):
        cooling_plate_section = CoolingPlatePart(thermal_resistor, pipe, fluid)
    elif issubclass(pipe.__class__, hy.pipes.JunctionPipe):
        cooling_plate_section = CoolingPlateJunction(pipe, fluid)
    return cooling_plate_section
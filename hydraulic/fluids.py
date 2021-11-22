#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fluids database
"""
from copy import copy

class Fluid:
    def __init__(self, rho, nu, heat_capacity, name=''):
        self.rho = rho
        self.nu = nu
        self.heat_capacity = heat_capacity
        self.name = name

    def Dict(self):
        return copy(self.__dict__)

    @classmethod
    def dict_to_object(cls, dict_):
        fluid = cls(dict_['rhos'], dict_['nu'],
                    dict_['heat_capacity'], dict_['name'])
        return fluid
        
water = Fluid(1000, 0.0010518, 4185, 'Water')

# Approximative figures
sae10 = Fluid(900, 0.1, 2000, 'SAE10 oil')
sae20 = Fluid(900, 0.2, 2000, 'SAE20 oil')
sae30 = Fluid(900, 0.5, 2000, 'SAE30 oil')
sae40 = Fluid(900, 0.8, 2000, 'SAE40 oil')
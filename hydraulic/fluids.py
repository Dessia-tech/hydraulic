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
    def DictToObject(cls, dict_):
        fluid = cls(dict_['rhos'], dict_['nu'],
                    dict_['heat_capacity'], dict_['name'])
        return fluid
        
water = Fluid(1000, 0.0010518, 4185, 'Water')

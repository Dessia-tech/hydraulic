#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fluids database
"""

class Fluid:
    def __init__(self, rho, nu, name=''):
        self.rho = rho
        self.nu = nu
        self.name = name
        
water = Fluid(1000, 0.0010518, 'Water')

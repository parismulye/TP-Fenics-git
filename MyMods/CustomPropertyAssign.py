#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 17:51:04 2020

@author: simone
"""
import sys
from   dolfin import FunctionSpace, Function

def AssignSubPropertiesLame(mesh, subdomains, E, nu):
    # To assign different coefficients to subdomains, we create a function space of order 0 (coefficients):
    V0 = FunctionSpace(mesh, 'DG', 0)
    lmbda = Function(V0)
    mu = Function(V0)  
    lmbda_array = (E*nu)/((1+nu)*(1-2*nu))
    mu_array = E/(2*(1+nu))    
    try:
        for cell_no in range(len(subdomains)):
            subdomain_no = subdomains[cell_no]
            lmbda.vector()[cell_no] = lmbda_array[subdomain_no]
            mu.vector()[cell_no] = mu_array[subdomain_no]
    except IndexError:
        print('Make sure you have subdomains and k values in equal number')
        sys.exit()        
    return (lmbda, mu)

def AssignSubPropertiesScalar(mesh, subdomains, k_values):
    V0 = FunctionSpace(mesh, 'DG', 0)
    k = Function(V0)
    try:
        for cell_no in range(len(subdomains)):
            subdomain_no = subdomains[cell_no]
            k.vector()[cell_no] = k_values[subdomain_no]
    except IndexError:
        print('Make sure you have subdomains and k values in equal number')
        sys.exit()
    return k
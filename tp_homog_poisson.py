#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 17:51:04 2020

@author: simone
"""

#  POISSON EQUATION

import numpy as np
import matplotlib.pyplot as plt
from   dolfin import *
from   MyMods import FilePrep as fpr
from   MyMods import CustomPropertyAssign as cpa


fname       = fpr.OpenGMSHGeom()
mesh        = Mesh(fname + ".xml")
subdomains  = MeshFunction("size_t", mesh, fname + "_physical_region.xml")      
boundaries  = MeshFunction("size_t", mesh, fname + "_facet_region.xml")

# Compute integration measures for this mesh:
ds          = Measure('ds', domain=mesh, subdomain_data=boundaries)
dx          = Measure('dx', domain=mesh, subdomain_data=subdomains)
volume      = 1

# To assign different coefficients k to subdomains:
k_values    = [0.1, 0.001]
k           = cpa.AssignSubPropertiesScalar(mesh, subdomains, k_values)


# u and v are scalar valued. Hence we create a suitable function space
V           = FunctionSpace(mesh, 'P', 1)
u           = TrialFunction(V)
v           = TestFunction(V)

# The left side of the weak form:
a = k*dot(grad(u), grad(v))*dx

# The left side of the weak form with source term:
f = Constant(0.0)
L = f*v*dx


# BCs scheme (in the provided geometries, not always! Depends on the id assigned in GMSH):
#        
#      _ _4_ _
#     |       |
#    1|       |3
#     |_ _ _ _|
#         2
# 

bcs = [
       DirichletBC(V, 0.0, boundaries, 1),                                      
       DirichletBC(V, 1.0, boundaries, 3)
       ]

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# ------------------------------ POST PROCESSING ------------------------------

W = VectorFunctionSpace(mesh, 'P', 2)
Q = project(k*grad(u), W)
qx, qy = Q.split(deepcopy=True)  
qx_nodal_values = np.array(qx.vector())
qy_nodal_values = np.array(qy.vector())

imposed_grad_u = 1
k_hom_x = np.mean(qx_nodal_values) / imposed_grad_u
k_hom_y = np.mean(qy_nodal_values) / imposed_grad_u

print('The apparent value of kx is: {}'.format(k_hom_x))
print('The apparent value of ky is: {}'.format(k_hom_y))

# ------------------------------- PLOTTING ------------------------------------

plt.figure(1)
plot(mesh, linewidth=0.5, color='black')
p = plot(subdomains, title='Subdomains')
plt.colorbar(p)

plt.figure(2)
p = plot(u, title='Solution U')
plt.colorbar(p)

plt.figure(3)
plot(subdomains)
p = plot(Q, title='Flux')
plt.colorbar(p)

fpr.ExportToParaview(fname, mesh, u, Flux=Q)
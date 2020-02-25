#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 17:51:04 2020

@author: simone
"""

# LINEAR MECHANICAL PROBLEM

import numpy as np
import matplotlib.pyplot as plt
from   dolfin import *
from   MyMods import FilePrep as fpr
from   MyMods import CustomBC as cbc
from   MyMods import CustomPropertyAssign as cpa


fname       = fpr.OpenGMSHGeom()
mesh        = Mesh(fname + ".xml")
subdomains  = MeshFunction("size_t", mesh, fname + "_physical_region.xml")      
boundaries  = MeshFunction("size_t", mesh, fname + "_facet_region.xml")

# Compute integration measures for this mesh:
ds          = Measure('ds', domain=mesh, subdomain_data=boundaries)
dx          = Measure('dx', domain=mesh, subdomain_data=subdomains)
volume      = 1

# Material parameters definition: we call a custom function that will create discontinuous functions (lmbda, mu) defined on subdomains
E           = np.array([10000, 10000])                                         
nu          = np.array([0.3, 0.3])
lmbda, mu   = cpa.AssignSubPropertiesLame(mesh, subdomains, E, nu)

# u and v are vector valued. Hence we create a suitable function space
Ve          = VectorElement("CG", mesh.ufl_cell(), 2)
V           = FunctionSpace(mesh, Ve)
u           = TrialFunction(V)
v           = TestFunction(V)

# Define strain and Hooke's constitutive law:
def epsilon(u):
    return sym(grad(u))

def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(2) + 2*mu*epsilon(u)

# The left side of the weak form:
a = inner(sigma(u), epsilon(v))*dx

# edge Neumann condition: uniform traction on boundary denoted by ds(ID). In this example, edge 3.
t = Constant((0.0, 0.0))                                                       
l_bound = dot(t, v)*ds(3)

# body forces:
f = Constant([0.0,0.0])
l_body = dot(f, v)*dx

# The right side of the weak form:
L = l_body + l_bound

# BCs scheme (in the provided geometries, not always! Depends on the id assigned in GMSH):
#        
#      _ _4_ _
#     |       |
#    1|       |3
#     |_ _ _ _|
#         2
# 

#11:
bcs = [
        DirichletBC(V, Constant((0,0)), cbc.bot_left_vertex, method='pointwise'),
        DirichletBC(V.sub(1), Constant(0), boundaries, 2),
        DirichletBC(V.sub(1), Constant(0.1), boundaries, 4)
      ]

# Compute solution:
u = Function(V)
solve(a == L, u, bcs)

# ------------------------------ POST PROCESSING ------------------------------

def Tensor2Voigt(s):
    return as_vector([s[0,0], s[1,1], s[0,1]])

sigma_voigt = Tensor2Voigt(sigma(u))
eps_voigt   = Tensor2Voigt(epsilon(u))

Sigma_avg   = np.zeros((3,))
Eps_avg     = np.zeros((3,))
for k in range(3):
    Sigma_avg[k]    = assemble(sum([sigma_voigt[k]*dx]))/volume
    Eps_avg[k]      = assemble(sum([eps_voigt[k]*dx]))/volume
    
print('Average stress in Voigt notation:', Sigma_avg)
print('Average strain in Voigt notation:', Eps_avg)

# ------------------------------- PLOTTING ------------------------------------

plt.figure(1)
plot(mesh, linewidth=0.5, color='black')
p = plot(subdomains, title='Subdomains')
plt.colorbar(p)

plt.figure(2)
p = plot(u, mode='displacement', title='Displacement U')
plt.colorbar(p)

# -------------------------- EXPORTING RESULTS --------------------------------

fpr.ExportToParaview(fname, mesh, u, Stress=sigma_voigt, Strain=eps_voigt)
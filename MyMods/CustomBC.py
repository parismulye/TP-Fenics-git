#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dolfin import SubDomain, DOLFIN_EPS



class Periodic_TB(SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[1] < DOLFIN_EPS and x[1] > -DOLFIN_EPS and on_boundary)
    def map(self, x, y):
        y[1] = x[1] - 1.0
        y[0] = x[0]

class Periodic_LR(SubDomain):
    def inside(self, x, on_boundary):
        return bool(x[0] < DOLFIN_EPS and x[0] > -DOLFIN_EPS and on_boundary)
    def map(self, x, y):
        y[0] = x[0] - 1.
        y[1] = x[1]
        

# define the 0,0 point:
def bot_left_vertex(x,on_boundary):
    return (abs(x[0]) < DOLFIN_EPS) and (abs(x[1]) < DOLFIN_EPS)

# define the 1,0 point:
def bot_right_vertex(x,on_boundary):
    return (abs(x[0] - 1.) < DOLFIN_EPS) and (abs(x[1]) < DOLFIN_EPS)

# define the 1,1 point:
def top_right_vertex(x,on_boundary):
    return (abs(x[0] - 1.) < DOLFIN_EPS) and (abs(x[1] - 1.) < DOLFIN_EPS)

# define the 0,01 point:
def top_left_vertex(x,on_boundary):
    return (abs(x[0]) < DOLFIN_EPS) and (abs(x[1] - 1.) < DOLFIN_EPS)
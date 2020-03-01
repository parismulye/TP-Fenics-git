#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess

from tkinter import Tk, messagebox
from tkinter.filedialog import askopenfilename
from dolfin import XDMFFile, VectorFunctionSpace, Function, project


def OpenGMSHGeom():
    Tk().withdraw()
    messagebox.showinfo("Information","Select the GMSH geometry file (.geo)")
    file_path = askopenfilename(filetypes =[('GMSH Files', '*.geo')])           # show an "Open" dialog box and return the path to the selected file
    fname = os.path.splitext(os.path.basename(file_path))[0]
    path = os.getcwd()
    print ('The current directory is %s' % path)
    dir_path = './mesh_files/{}'.format(fname)
    try:
        os.makedirs(dir_path)
    except OSError:
        print ("Creation of the mesh directory %s failed. Maybe it already exists..." % dir_path)
    else:
        print ("Successfully created the mesh directory %s " % dir_path)
    os.system('gmsh {}.geo -format msh2 -2 -o mesh_files/{}/{}.msh -v 3'.format(fname, fname, fname))
    os.chdir('mesh_files/{}/'.format(fname))
    subprocess.run(['dolfin-convert', '{}.msh'.format(fname), '{}.xml'.format(fname)], stdout=subprocess.DEVNULL)
    return fname



def ExportToParaview(fname, mesh, u, **kwargs):
    file_results = XDMFFile('{}_results.xdmf'.format(fname))
    file_results.parameters["flush_output"] = True
    file_results.parameters["functions_share_mesh"] = True
    V2 = VectorFunctionSpace(mesh, "CG", degree=1, dim=2)
    V3 = VectorFunctionSpace(mesh, "CG", degree=1, dim=3)
    for key,val in kwargs.items():
        if len(val) == 3:
            res = Function(V3, name=key)
        elif len(val) == 2:
            res = Function(V2, name=key)
        res.assign(project( val ) )
        file_results.write(res, 0.)
    u.rename('u', 'displacement')
    file_results.write(u, 0.)

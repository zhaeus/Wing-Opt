# -*- coding: utf-8 -*-
"""
Script for processing and graphing the geometry for the wing optimisation script
"""
import math 
import numpy as np 
import matplotlib.pyplot as plt 
from openpyxl import Workbook
import os

## Geometric parameters 

span = 9.108 #m, wingspan
semispan = span/2
planform_area = 15.2212 #m2
chord = planform_area/span #

## Aerofoil 

naca = np.loadtxt("naca0313.txt")
naca_len = len(naca)
half_len = int(naca_len/2) + 1
naca *= chord # Scaling loop from unit length to average chord length

left_right_OML_x = naca[:,0]
left_right_OML_z = naca[:,1]

for column in 0, 1: # Reversing direction of points in text input
    naca[half_len:,column] = naca[half_len:,column][::-1] 
    
# Defining inner and outer skin surfaces
x_outer = naca[:,0]
x_outer = np.append(x_outer,x_outer[0])
z_outer = naca[:,1]
z_outer = np.append(z_outer,z_outer[0])

## For use in Excel 

# for column in 0,1: 
#     print('start')
#     for row in range(naca_len):
#         print(f"{naca[row,column]}")



x_outer_upper = x_outer[:half_len]
x_outer_lower = x_outer[half_len:]
z_outer_upper = z_outer[:half_len]
z_outer_lower = z_outer[half_len:]



# xsection_area = np.abs(green_area(x_outer,z_outer))

# def inset_offset(xvec,yvec,thickness):
#     lenx = len(xvec)
#     leny = len(yvec)
#     if lenx != leny:
#         raise RuntimeError('Vectors must have same length')
#     if lenx < 3:
#         raise RuntimeError('Vectors must form closed polygon')
#     if xvec[0] != xvec[-1] or yvec[0] != yvec[-1]:
#         np.append(xvec,xvec[0])
#         np.append(yvec,yvec[0])
#     xvec_set = np.zeros(lenx+1)
#     yvec_set = np.zeros(leny+1)
#     for i in range(1,lenx-1):
#         try:
#             m = (yvec[i+1]-yvec[i-1])/(xvec[i+1]-xvec[i-1])
#             # m = 0.5*((yvec[i]-yvec[i-1])*(xvec[i+1]-xvec[i]) + (yvec[i+1]-yvec[i])*(xvec[i]-xvec[i-1]))
#             mdash = -1/m
#         except:
#             mdash = 0
#         sign = (2*bool(yvec[i]>yvec[i-1])-1)
#         if mdash == 0:
#             xvec_set[i] = sign*thickness + xvec[i]
#             yvec_set[i] = yvec[i]
#         xvec_set[i] = np.sqrt((thickness**2)/((mdash**2)+1))*sign + xvec[i]
#         yvec_set[i] = (xvec_set[i] - xvec[i])*mdash + yvec[i] 
#     set_curve = np.vstack((xvec_set.T,yvec_set.T))
#     return set_curve

# def plot_aerofoil():
#     _, ax = plt.subplots()
#     ax.plot(naca[:,0], naca[:,1])
#     plt.title(f"Area = {xsection_area:.3f} m2")
#     plt.ylim(-1, 1)
#     plt.show()
    
# skin_thickness = 0.0015 #m
# inner_surface = inset_offset(x_outer,z_outer,skin_thickness)
# x_inner = inner_surface[0,:]
# z_inner = inner_surface[1,:]

# def plot_innerouter_surface():
#     _, ax = plt.subplots()
#     ax.plot(x_inner,z_inner,color='orange',label='Inner')
#     ax.plot(x_outer,z_outer,color='blue',label='Outer - NACA0313')
#     plt.title('Inner surface and outer surface')
#     # plt.ylim(-chord/2, chord/2)
#     plt.ylim(-0.2,0.2)
#     plt.legend()
#     plt.show()

# # Method 1: Surface
# def plot_aerofoil_surface():
#     fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    
#     x = x_outer_upper 
    
#     span_num = 10
#     y = np.linspace(0,semispan,span_num)
    
#     xx, yy = np.meshgrid(x,y)
    
#     z = z_outer_upper
#     zz,_ = np.meshgrid(z,y)
    
#     ax.plot_surface(xx, yy, zz,color='red')
    
#     ax.set_title('NACA0313 Airfoil')
#     ax.set_zlim(-semispan/3,semispan/3)
#     ax.set_xlim(-semispan,semispan)
#     plt.show()

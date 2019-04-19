#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 12:01:54 2018

@author: abauville

MIT License

Copyright (c) 2019 abauville

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from FaultPlot import plotFaultArrow
import numpy as np
from numpy import array as arr
from numpy import sin,cos,pi,tan
import matplotlib.pyplot as plt

def intersection(segment1_x, segment1_y, segment2_x, segment2_y):
    p = arr([segment1_x[0],segment1_y[0]])
    r = arr([segment1_x[1]-segment1_x[0],segment1_y[1]-segment1_y[0]])
    
    q = arr([segment2_x[0],segment2_y[0]])
    s = arr([segment2_x[1]-segment2_x[0],segment2_y[1]-segment2_y[0]])

    t = np.cross( (q-p) , (s/np.cross(r,s)) ) 
    

    return [(p+t*r),t]

def plotArrow(x,y,theta=0.0,length=1.0,bodyWidth=0.0015,headLength=0.25,headWidth=0.25/3.0,color="k",style='single'):
    
    #        sl = 0.15
    
    
    try: # x and y have length 2
        theta=np.arctan2((y[1]-y[0]),(x[1]-x[0]))
#        if []
#        print(str(theta*180.0/pi))
#        if np.abs(theta)>1.0*np.pi:
#            theta-=1.0*np.pi
        L = np.sqrt((y[1]-y[0])**2 + (x[1]-x[0])**2)
        x = x[1]
        y = y[1]
    except IndexError: # x and y or given as scalars
        OK = 1
        L = length
    
    aHL = headLength# = L/4.0 # arrow head length
    aHW = headWidth# L/12.0
    aBL = L-aHL
    aBW = bodyWidth#aHW/3.0
    
    a_x = np.array([-aHL, -aHL, -aHL-aBL, -aHL-aBL, -aHL, -aHL, 0.0]) # arrow points x
    a_y = np.array([+aHW, +aBW,   +aBW  ,   -aBW  , -aBW, -aHW, 0.0]) # arrow points x
    
    
    
    
    
    
    
    if style=='single':
        rot_a_x = np.cos(theta)*a_x - np.sin(theta)*a_y
        rot_a_y = np.sin(theta)*a_x + np.cos(theta)*a_y
        plt.fill(x+rot_a_x,y+rot_a_y,lineStyle='None',color=color)
    elif style=='opposing':
        rot_a_x = np.cos(theta)*a_x - np.sin(theta)*a_y
        rot_a_y = np.sin(theta)*a_x + np.cos(theta)*a_y
        plt.fill(x+rot_a_x,y+rot_a_y,lineStyle='None',color=color)
        plt.fill(x-rot_a_x,y-rot_a_y,lineStyle='None',color=color)
    elif style=='double':
        
        a_x+=L
        a_x/=2.0
        
        rot_a_x = np.cos(theta)*a_x - np.sin(theta)*a_y
        rot_a_y = np.sin(theta)*a_x + np.cos(theta)*a_y
        x += np.cos(theta)*L/2.0 - np.sin(theta)*0
        y += np.cos(theta)*0 - np.sin(theta)*L/2.0
        plt.fill(x+rot_a_x,y+rot_a_y,lineStyle='None',color=color)
        plt.fill(x-rot_a_x,y-rot_a_y,lineStyle='None',color=color) 
    else:
        raise ValueError ("Unknwon style. style can only take the values 'single', 'opposing', 'double'.")

def plotWedge(taper,enveloppe="lower",beta=0.0,
              origin=arr([0.0,0.0]),
              fx0_list_a = arr([.2, .4, .6, .8]),
              fy0_list_a = arr([0.0]),
              fx0_list_b = arr([.2, .4, .6, .8]),
              fy0_list_b = arr([0.0]),
              faultPos = 0.33,
              plotFaults=True, plotWedge=True,plotFaultsArrow=True,plotStress=True,
              sx0=0.9,sy0=0.1,
              colorWedge=[.9,.9,.95,0.5],colorSurface="k",colorBase="k",
              sl=0.075,s_arrowHeadLength = 1.0/2.0, s_arrowHeadWidth = 1.0/6.0, s_arrowBodyWidth = 1.0/3.0,
              f_arrowLength=.02, f_arrowHeadLength = .5,f_arrowSpacing=0.01,
              lineWidthSurface=1.0, lineWidthBase=1.0,
              colorFaults='r',lineWidthFaults=.5,
              back_x = 1.0):

    
    
    if enveloppe=='lower':
        alpha = taper.findAlpha(beta,"lower")
        psi = taper.psi_bmin
    elif enveloppe=='upper':
        alpha = taper.findAlpha(beta,"upper")
        psi = taper.psi_bmax
    elif enveloppe=='average':
        alpha = taper.findAlpha(beta,"upper")
        alpha += taper.findAlpha(beta,"lower")
        alpha /= 2.0
        psi = (taper.psi_bmax+taper.psi_bmin)/2.0
    else:
        raise ValueError("Enveloppe should be 'upper' or 'lower'")
        
    tAngle = alpha+beta
#    print(tAngle)
        
    # Plot wedge outline
    # ===================================
    surf_x = origin[0] + arr([0.0,1.0])
    surf_y = origin[1] + arr([0.0, sin(tAngle)*(2.0-cos(tAngle))])
    back_x = origin[0] + arr([back_x,back_x])
    back_y = origin[1] + arr([0.0,sin(tAngle)*(2.0-cos(tAngle))])
    base_x = origin[0] + arr([0.0,1.0])
    base_y = origin[1] + arr([0.0,0.0])
    if plotWedge:
        plt.fill(np.concatenate((surf_x, arr([origin[0]+1.0]))),np.concatenate((surf_y, arr([origin[1]]))),color=colorWedge)
        plt.plot([origin[0], origin[0]+1.0], [origin[1],origin[1]],color=colorBase,lineWidth=lineWidthBase)
        plt.plot(surf_x,surf_y,color=colorSurface,lineWidth=lineWidthSurface)
    
    
        # Plot stress orientation
        # ===================================
    if plotStress: 
        sx0 = origin[0] + sx0
        sy0 = origin[1] + sy0

        aHL = s_arrowHeadLength*sl
        aHW = s_arrowHeadWidth*sl
        aBW = aHW*s_arrowBodyWidth

        plotArrow(sx0,sy0,psi,length=sl,headLength=aHL,headWidth=aHW, bodyWidth=aBW,style='opposing')
    if plotFaults:
        # Plot faults
        # ===================================
        
        ca = 30.0*pi/180.0 # Coulomb angle
        
        fa_list = psi + np.array([-ca, +ca])
        #fa_list = psi + np.array([+ca])
        
        
        

        for fa in fa_list:
            if fa==fa_list[0]:
                fx0_list = origin[0] + fx0_list_a
                fy0_list = origin[1] + fy0_list_a
            else:
                fx0_list = origin[0] + fx0_list_b
                fy0_list = origin[1] + fy0_list_b
                
            for fx0 in fx0_list:
                for fy0 in fy0_list:
                
                    # check intersection
                    fx = arr([fx0-cos(fa),fx0+cos(fa)])
                    fy = arr([fy0-sin(fa),fy0+sin(fa)])
                    
                    it, t = intersection(fx,fy, surf_x,surf_y)
                    if (fy[0]-origin[1]<(fx[0]-origin[0])*tan(tAngle)):
                        fx = arr([fx[0],it[0]])
                        fy = arr([fy[0],it[1]])
                    else:
                        fx = arr([it[0],fx[1]])
                        fy = arr([it[1],fy[1]])
                        
                    it, t = intersection(fx,fy, base_x,base_y)
                    if (t>=0.0 and t<=1.0):
                        if (fy[0]-origin[1]>0.0):
                            fx = arr([fx[0],it[0]])
                            fy = arr([fy[0],it[1]])
                        else:
                            fx = arr([it[0],fx[1]])
                            fy = arr([it[1],fy[1]])

                    it, t = intersection(fx,fy, back_x,back_y)
                    if (t>0.0 and t<1.0):
                        fx = arr([fx[0],it[0]])
                        fy = arr([fy[0],it[1]])                 
                        
                        
                
                    plt.plot( fx , fy , '-', color=colorFaults,linewidth=0.75)
                    
                    if plotFaultsArrow:
                        if fa>psi:
                            u = faultPos
                            sense=1
                        else:
                            u = 1.0-faultPos
                            sense=0
                        x = (1.0-u)*fx[0] + u*fx[1]
                        y = (1.0-u)*fy[0] + u*fy[1]
                        plotFaultArrow(x,y,fa, L=f_arrowLength, headL = f_arrowHeadLength, sense=sense, spacing=f_arrowSpacing, color=colorFaults,linewidth = lineWidthFaults)
                        plotFaultArrow(x,y,fa, L=f_arrowLength, headL = f_arrowHeadLength, sense=sense, spacing=-f_arrowSpacing, color=colorFaults,linewidth = lineWidthFaults)
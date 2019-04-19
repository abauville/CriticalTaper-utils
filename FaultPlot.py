#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 16:52:27 2017

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

import numpy as np
from numpy import sin, cos, tan, arcsin, arccos, arctan, pi
import matplotlib.pyplot as plt

degree = 180.0/pi



# =============================================================================
#
#                    Functions for fault diagram plotting
#
# =============================================================================
def plotFaultArrow(x,y,theta, L=1, sense=0, spacing=0.1, color="r",angleHead=20.0 * pi/180.0,headL = .5,ax=plt,linewidth = 1):
    # sense 0: sinistral, 1:dextral
    x = x + sin(theta)*spacing
    y = y - cos(theta)*spacing
    segment = np.array((-1,1)) * L
    segmentHead = np.array((0,2)) * L*headL
    
    ax.plot(x+cos(theta)*segment, y   + sin(theta)*segment,color=color,linewidth=linewidth)
    if ((spacing>0) & (sense == 0)):
        ax.plot(x+L*cos(theta) - cos(angleHead+theta)*(segmentHead), y+L*sin(theta) - sin(angleHead+theta)*(segmentHead),color=color,linewidth=linewidth)
    elif ((spacing>0) & (sense == 1)):
        ax.plot(x-L*cos(theta) + cos(-angleHead+theta)*(segmentHead), y-L*sin(theta) + sin(-angleHead+theta)*(segmentHead),color=color,linewidth=linewidth)
    elif ((spacing<=0) & (sense == 0)):
        ax.plot(x-L*cos(theta) + cos(angleHead+theta)*(segmentHead), y-L*sin(theta) + sin(angleHead+theta)*(segmentHead),color=color,linewidth=linewidth)
    elif ((spacing<=0) & (sense == 1)):
        ax.plot(x+L*cos(theta) - cos(-angleHead+theta)*(segmentHead), y+L*sin(theta) - sin(-angleHead+theta)*(segmentHead),color=color,linewidth=linewidth)
    else:
        raise ValueError("sense must be 0 or 1")
    
    
def plotArrow(x,y,theta, L=1, color="r", sense=0, angleHead=20.0 * pi/180.0,headL = .5,ax=plt,linewidth = 1):
    # sense 0: sinistral, 1:dextral
    x = x# + sin(theta)*spacing
    y = y# - cos(theta)*spacing
    segment = np.array((-1,1)) * L
    segmentHead = np.array((0,2)) * L*headL
    
    ax.plot(x+cos(theta)*segment, y   + sin(theta)*segment,color=color,linewidth=linewidth)

    if sense == 0:
        ax.plot(x+L*cos(theta) - cos(angleHead+theta)*(segmentHead), y+L*sin(theta) - sin(angleHead+theta)*(segmentHead),color=color,linewidth=linewidth)
#    elif ((spacing>0) & (sense == 1)):

        ax.plot(x+L*cos(theta) - cos(-angleHead+theta)*(segmentHead), y+L*sin(theta) - sin(-angleHead+theta)*(segmentHead),color=color,linewidth=linewidth)
    elif sense == 1:
#        ax.plot(x-L*cos(theta) + cos(-angleHead+theta)*(segmentHead), y-L*sin(theta) + sin(-angleHead+theta)*(segmentHead),color=color,linewidth=linewidth)
        ax.plot(x-L*cos(theta) + cos(angleHead+theta)*(segmentHead), y-L*sin(theta) + sin(angleHead+theta)*(segmentHead),color=color,linewidth=linewidth)
#    elif ((spacing<=0) & (sense == 1)):
        ax.plot(x-L*cos(theta) + cos(-angleHead+theta)*(segmentHead), y-L*sin(theta) + sin(-angleHead+theta)*(segmentHead),color=color,linewidth=linewidth)
    else:
        raise ValueError("sense must be 0 or 1")
    
    

def plotFaultDiagram(Tau,psi, L=1,colorFault="r",colorSigma1="b",Larrow=.15,PosArrow=.66,angleHeadArrow=20.0 * pi/180.0, spacing=.1,ax=plt,refAxes=0,faultLinewidth=1,arrowLinewidth=1,sigma1Linewidth=1,polar=0):
    segment = np.array((-1,1)) * L
    thetaA = psi+30*pi/180
    thetaB = psi-30*pi/180
    
    if (refAxes==0):
        Tau = Tau
        psiPos = psi*degree
    else:
        Tau = ( (Tau - refAxes.axis()[0])/(refAxes.axis()[1]-refAxes.axis()[0]) - ax.axis()[0] ) * (ax.axis()[1]-ax.axis()[0])
        psiPos = ( (psi*degree - refAxes.axis()[2])/(refAxes.axis()[3]-refAxes.axis()[2]) - ax.axis()[2] ) * (ax.axis()[3]-ax.axis()[2])
        
    
    if (polar==0):
        # Sigma1 dir
        ax.plot(Tau+cos(psi)*segment, psiPos  + sin(psi)*segment,color=colorSigma1,linewidth=sigma1Linewidth)
        
        # Faults
        ax.plot(Tau+cos(thetaA)*segment, psiPos   + sin(thetaA)*segment,color=[.8,.5,.2],linewidth=faultLinewidth)
        ax.plot(Tau+cos(thetaB)*segment, psiPos   + sin(thetaB)*segment,color=[.6,.3,.6],linewidth=faultLinewidth)
    else:
        # Sigma1 dir
        ax.plot(Tau*cos(psi)+cos(psi)*segment, Tau*sin(psi)  + sin(psi)*segment,color=colorSigma1,linewidth=sigma1Linewidth)
        
        # Faults
        ax.plot(Tau*cos(psi)+cos(thetaA)*segment, Tau*sin(psi)   + sin(thetaA)*segment,color=[.8,.5,.2],linewidth=faultLinewidth)
        ax.plot(Tau*cos(psi)+cos(thetaB)*segment, Tau*sin(psi)   + sin(thetaB)*segment,color=[.6,.3,.6],linewidth=faultLinewidth)
        
    
#    ax.plot(Tau+cos(thetaA)*segment, psiPos   + sin(thetaA)*segment,color=colorFault,linewidth=faultLinewidth)
#    ax.plot(Tau+cos(thetaB)*segment, psiPos   + sin(thetaB)*segment,color=colorFault,linewidth=faultLinewidth)
    
    # Arrows

    # All arrows   
#    plotFaultArrow(Tau-cos(thetaA)*PosArrow*L,psiPos-sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing=spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#    plotFaultArrow(Tau-cos(thetaA)*PosArrow*L,psiPos-sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#    plotFaultArrow(Tau-cos(thetaB)*PosArrow*L,psiPos-sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#    plotFaultArrow(Tau-cos(thetaB)*PosArrow*L,psiPos-sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing= spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#    
#    plotFaultArrow(Tau+cos(thetaA)*PosArrow*L,psiPos+sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing= spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#    plotFaultArrow(Tau+cos(thetaA)*PosArrow*L,psiPos+sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#    plotFaultArrow(Tau+cos(thetaB)*PosArrow*L,psiPos+sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#    plotFaultArrow(Tau+cos(thetaB)*PosArrow*L,psiPos+sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing= spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#    
#
#    # Outer arrows only
#    plotFaultArrow(Tau-cos(thetaA)*PosArrow*L,psiPos-sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing=spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#    plotFaultArrow(Tau-cos(thetaB)*PosArrow*L,psiPos-sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#
#    plotFaultArrow(Tau+cos(thetaA)*PosArrow*L,psiPos+sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#    plotFaultArrow(Tau+cos(thetaB)*PosArrow*L,psiPos+sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing= spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
#    
    if (polar==0):
        # Inner arrows only
        plotFaultArrow(Tau-cos(thetaA)*PosArrow*L,psiPos-sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
        plotFaultArrow(Tau-cos(thetaB)*PosArrow*L,psiPos-sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing= spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
        
        plotFaultArrow(Tau+cos(thetaA)*PosArrow*L,psiPos+sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing= spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
        plotFaultArrow(Tau+cos(thetaB)*PosArrow*L,psiPos+sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
    
    else:
        # Inner arrows only
        plotFaultArrow(Tau*cos(psi)-cos(thetaA)*PosArrow*L,Tau*sin(psi)-sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
        plotFaultArrow(Tau*cos(psi)-cos(thetaB)*PosArrow*L,Tau*sin(psi)-sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing= spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
        
        plotFaultArrow(Tau*cos(psi)+cos(thetaA)*PosArrow*L,Tau*sin(psi)+sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing= spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
        plotFaultArrow(Tau*cos(psi)+cos(thetaB)*PosArrow*L,Tau*sin(psi)+sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)







def plotFaultDiagram2(Tau,psi,x,y,L=1,colorFault="r",colorSigma1="b",Larrow=.15,PosArrow=.66,angleHeadArrow=20.0 * pi/180.0, spacing=.1,ax=plt,refAxes=0,faultLinewidth=1,arrowLinewidth=1,sigma1Linewidth=1,polar=0):
    segment = np.array((-1,1)) * L
    thetaA = psi+30*pi/180
    thetaB = psi-30*pi/180
    
    
    if (refAxes==0):
        Tau = Tau
#        psiPos = psi*degree
    else:
        Tau = ( (Tau - refAxes.axis()[0])/(refAxes.axis()[1]-refAxes.axis()[0]) - ax.axis()[0] ) * (ax.axis()[1]-ax.axis()[0])
#        psiPos = ( (psi*degree - refAxes.axis()[2])/(refAxes.axis()[3]-refAxes.axis()[2]) - ax.axis()[2] ) * (ax.axis()[3]-ax.axis()[2])

    # Sigma1 dir
#    ax.plot(Tau*cos(psi)+cos(psi)*segment, Tau*sin(psi)  + sin(psi)*segment,color=colorSigma1,linewidth=sigma1Linewidth)
    ax.plot(x+cos(psi)*segment, y+ sin(psi)*segment,':',color=colorSigma1,linewidth=sigma1Linewidth)
    
    # Faults
    ax.plot(x+cos(thetaA)*segment, y  + sin(thetaA)*segment,color=[.8,.5,.2],linewidth=faultLinewidth)
    ax.plot(x+cos(thetaB)*segment, y  + sin(thetaB)*segment,color=[.6,.3,.6],linewidth=faultLinewidth)
    

    # Inner arrows only
    plotFaultArrow(x-cos(thetaA)*PosArrow*L,y-sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
    plotFaultArrow(x-cos(thetaB)*PosArrow*L,y-sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing= spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
    
    plotFaultArrow(x+cos(thetaA)*PosArrow*L,y+sin(thetaA)*PosArrow*L,thetaA,sense=1,spacing= spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)
    plotFaultArrow(x+cos(thetaB)*PosArrow*L,y+sin(thetaB)*PosArrow*L,thetaB,sense=0,spacing=-spacing,L=Larrow,ax=ax,linewidth=arrowLinewidth,angleHead=angleHeadArrow)















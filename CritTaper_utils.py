#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 20:49:33 2018

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

# Critical Taper utilities
import numpy as np
from numpy import sin, cos, tan, pi, arctan, arcsin
import matplotlib.pyplot as plt

class Taper():
        
    def __init__(self, phi=30.0*np.pi/180.0, phi_b=30.0*np.pi/180.0-1e-6,
                 Lambda=0.0, Lambda_b=0.0,
                 rho=2500.0, rho_w=1000.0):
        self.phi = phi
        self.phi_b = phi_b
        self.Lambda = Lambda
        self.Lambda_b = Lambda_b
        
        self.rho = rho
        self.rho_w = rho_w
        
        self.alpha_all = []
        self.beta_all = []
        self.psi_bmin = 0
        self.psi_bmax = 0


    
    
    
    # =========================================================================        
    #               Function to compute the graph alpha vs beta              
    def computeAlphaVsBeta(self, step0=0.01, n="default",alphaMin="default",alphaMax="default"):#
        # This functions constructs the critical taper envelope (CTE) using the 
        # analytical solution of Lehner (1986).
        # The nomenclature of angles follows Lehner (1986).
        # The solution has been benchmark against the graphical method of Dahlen (1984).
        #
        # If you use this function please cite this repository and the original paper:
        # cite 1
        # cite 2
        #
        # References:
        # Dahlen, F. A. (1984). Noncohesive critical Coulomb wedges: An exact solution. 
        #   Journal of Geophysical Research: Solid Earth, 89(B12):10125–10133.
        # 
        # Lehner, F. K. (1986). Comments on "Noncohesive critical Coulomb wedges: 
        #   an exact solution" by F. A Dahlen. Journal of Geophysical Research, 
        #   91(B1):793-796
        # 
        # Returns a 2*N matrix with [alpha;beta]; N=2*n
        # unit is radian
        # phi = internal friction angle
        # mu_b = basal coefficient of friction = ()
        # alpha = surface angle
        # beta = base angle
        
        # psi   = stress orientation
        # psi_b = stress orientation with respect to the base direction
        # psi_0 = stress orientation with respect to the surface direction
        
         
        phi = self.phi
        phi_b = self.phi_b
        Lambda = self.Lambda
        Lambda_b = self.Lambda_b
        
        rho_w = self.rho_w
        rho = self.rho
        
        if alphaMin!="default":
            print("warning: the parameter alphaMin is deprecated. It will be ignored.")
        
        if alphaMax!="default":
            print("warning: the parameter alphaMax is deprecated. It will be ignored.")

        if n!="default":
            print("warning: the parameter n is deprecated. It will be ignored. Use step0 to specify the the resolution")
        Lambda_hydro = rho_w/rho
        
        Lambda_ov   = 1.0 - (1.0-Lambda  )/(1.0-Lambda_hydro)
        Lambda_b_ov = 1.0 - (1.0-Lambda_b)/(1.0-Lambda_hydro)
        
        
        mu = tan(phi)
        mu_b = tan(phi_b)
        mu_b_p = mu_b*((1.0-Lambda_b)/(1.0-Lambda))
        phi_b_p = arctan(mu_b_p)
        
        alpha_m = np.arctan( (1.0-Lambda_b_ov)*np.tan(phi_b))
#        beta_m  = -alpha_m
#        beta_c  = np.pi/4.0-phi_b_p/2.0
        
        alpha_max = np.arctan((1.0-Lambda_ov)*mu)        
        
        alpha_list = [alpha_m,alpha_max-1e-10,-alpha_m,-alpha_max+1e-10,alpha_m]
        
        
        
        alpha_all = np.array([])
        beta_all = np.array([])
       
        
        for i in range(len(alpha_list)-1):
        #for i in range(1):
            if alpha_list[i]<alpha_list[i+1]:
                step = np.abs(step0)
            else:
                step = -np.abs(step0)
            
            if np.abs(1.0-np.abs(alpha_list[i]/alpha_max))<1e-6: 
                alpha1 = alpha_list[i]
                alpha0 = alpha_list[i+1]
            else:
                alpha0 = alpha_list[i]
                alpha1 = alpha_list[i+1]
                    
            stepRef = step
            
            alpha = np.array([])
            new_alpha = alpha_list[i]

            # Create the vector alpha using a smaller around alpha_max than around alpha_m
            while np.abs(new_alpha-alpha_list[i])<np.abs(alpha_list[i+1]-alpha_list[i]):
                alpha = np.append(alpha,new_alpha)
                if alpha0 == alpha_list[i]:
                    step = stepRef*(1.1-np.abs((new_alpha+step-alpha0)/(alpha1-alpha0)))
                else:
                    step = stepRef*(1.1-np.abs((new_alpha-alpha0)/(alpha1-alpha0)))
                new_alpha = new_alpha + step


            # Compute beta
            alpha_p = arctan( 1.0/(1.0-Lambda_ov)*tan(alpha) )      
            
            theta = arcsin(sin(phi_b_p)/sin(phi))               # auxiliary friction angle
            gamma = arcsin(sin(alpha_p)/sin(phi))               # auxiliary surface  angle
                
            if i%2==0:   psi_0 = 0.5*(-gamma - alpha_p + pi)    # psi_0 = psi_0_a
            else:        psi_0 = 0.5*(+gamma - alpha_p)         # psi_0 = psi_0_p
            
            if i<2:      psi_b = 0.5*(-theta - phi_b_p + pi)    # psi_0 = pi/2-psi_b_1
            else:        psi_b = 0.5*(+theta - phi_b_p)         # psi_0 = psi_b_3
        
            beta =  psi_b-psi_0-alpha
            beta[beta<-pi/4.0] += pi
            print(beta)
            
            if (i==0):
                self.psi_bmax = psi_b
            elif (i==3):
                self.psi_bmin = psi_b
            
            # append the segment results to the total vector
            alpha_all = np.append(alpha_all,alpha)
            beta_all = np.append(beta_all,beta)
            
        # end segment loop
        self.beta_all = beta_all
        self.alpha_all = alpha_all


    #               Function to compute the graph alpha vs beta              
    # =========================================================================
    
    
    
    # =========================================================================    
    #                      Find alpha for a given beta                       
    def findAlpha(self, beta , enveloppe, tol=1e-2):
        # enveloppe can be "upper" or "lower"
        # Find alpha for a given beta        
#        if beta>pi:
#            print("Warning: beta should be in radians")
        
        IBetaMin = np.argmin(self.beta_all)
        IBetaMax = np.argmax(self.beta_all)
    
        I0 = min(IBetaMin, IBetaMax)
        I1 = max(IBetaMin, IBetaMax)
        beta_env1 = self.beta_all[I0:I1+1]
        beta_env2 = np.concatenate((self.beta_all[I1::],self.beta_all[0:I0+1]))
    
        alpha_env1 = self.alpha_all[I0:I1+1]
        alpha_env2 = np.concatenate((self.alpha_all[I1::],self.alpha_all[0:I0+1]))
        Ienv1 = np.argmin(np.abs(beta_env1-beta))
        Ienv2 = np.argmin(np.abs(beta_env2-beta))
        
        if enveloppe=="upper":
            alpha = np.max([alpha_env1[Ienv1],alpha_env2[Ienv2]]) # min or max to get the lower or upper enveloppe of stability
        elif enveloppe=="lower":
            alpha = np.min([alpha_env1[Ienv1],alpha_env2[Ienv2]]) # min or max to get the lower or upper enveloppe of stability
        elif enveloppe=="average":
            alpha  = ( alpha_env1[Ienv1]+alpha_env2[Ienv2] )/2.0# min or max to get the lower or upper enveloppe of stability            
        else:
            raise ValueError("unknown enveloppe '%s'." % enveloppe)

        return alpha
    #                      Find alpha for a given beta                          
    # =========================================================================
    
    
    
#            
#    # =========================================================================        
#    #               Function to compute the graph alpha vs beta             
##    Deprecated function kept for potential debugging of the new function
#    def computeAlphaVsBeta_Numerical(self, n=1001,alphaMin="default",alphaMax="default"):
        # Some trigo functions
#        def csc(x):
#            return 1.0/sin(x)
#        def sec(x):
#            return 1.0/cos(x)
#        def cot(x):
#            return 1.0/tan(x)
#        
#        # F function from Dahlen, 1984, eq (21)
#        def F_function(psi,phi):
#            return tan(2.0*psi) / (csc(phi)*sec(2.0*psi) - 1.0)
#        # Computes the graph of surface vs base angle based on the graphical method described in
#        # Dahlen, 1984, Noncohesive Critical Coulomb Wedges (paragraph Stable vs Unstable Wedges p. 10,128).
#        # 
#        # Returns a 2*N matrix with [alpha;beta]; N=2*n
#        # unit is radian
#        # phi = internal friction angle
#        # mu_b = basal coefficient of friction = ()
#        # alpha = surface angle
#        # beta = base angle
#        
#        # psi   = stress orientation
#        # psi_b = stress orientation with respect to the base direction
#        # psi_0 = stress orientation with respect to the surface direction
#        
#        phi = self.phi
#        phi_b = self.phi_b
#        Lambda = self.Lambda
#        Lambda_b = self.Lambda_b
#        
#        rho_w = self.rho_w
#        rho = self.rho
#        
#        if alphaMin=="default":
#            alphaMin = -phi+1e-4
##        else: alphaMin = alphaMin
#            
#        if alphaMax=="default":
#            alphaMax = +phi-1e-4
#        
#        alphaP_list = np.linspace(alphaMin,alphaMax,n) #  note: alphaMax = phi
#        if np.where(alphaP_list==0)[0].shape[0] > 0: # if value 0 is found
#            # then increment n by one and reconstruct alphaP_list to avoid 0 (there is no solution for the value 0)
#            n=n+1
#            alphaP_list = np.linspace(-phi+1e-6,phi-1e-6,n) #  note: alphaMax = phi
#        
#        
#        
#        mu_bP = tan(phi_b) * (1.0-Lambda_b)/(1.0-Lambda)
#        if Lambda>0.0:
#            alpha_list = arctan(tan(alphaP_list) * (1.0-Lambda)/(1.0-rho_w/rho))
#    #        alpha_list = arctan(tan(alphaP_list) * (1.0-Lambda))
#        else:
#            alpha_list = alphaP_list
#                    
#        
#        psi = np.linspace(-pi/2.0,pi/2.0,n) 
#        psi_0    = np.zeros(n)
#        psi_02    = np.zeros(n)
#    
#        beta_Left = np.zeros(n)
#        beta_Right = np.zeros(n)
#        
#        F = F_function(psi,phi)
#    
#    
#        
#    
#        for i in range(n):
#            alphaP = alphaP_list[i]
#            alpha  = alpha_list[i]
#            # Find phi_b
#            I = np.where((F[:-1]-mu_bP)/(F[1:]-mu_bP) < 0.0 )[0]            
#            psi_bmin = (psi[I][0]+psi[I+1][0])/2.0
#            psi_bmax = (psi[I][1]+psi[I+1][1])/2.0
#            # Find phi_0
#            I = np.where((F[:-1]-tan(alphaP))/(F[1:]-tan(alphaP)) < 0.0 )[0][0]        
#            psi_0 = (psi[I]+psi[I+1])/2.0
#            I = np.where((F[:-1]-tan(alphaP))/(F[1:]-tan(alphaP)) < 0.0 )[0][1]
#            psi_02 = (psi[I]+psi[I+1])/2.0
#    
#    
#            ApBmin = psi_bmin-psi_0 # alpha + beta = critical taper
#            ApBmin2 = psi_bmin-psi_02 # alpha + beta = critical taper
#            ApBmax = psi_bmax-psi_0 # alpha + beta = critical taper
#            ApBmax2 = psi_bmax-psi_02 # alpha + beta = critical taper
#            
#            betaMin = ApBmin-alpha
#            betaMin2 = ApBmin2-alpha
#            betaMax = ApBmax-alpha
#            betaMax2 = ApBmax2-alpha
#        
#            if alphaP>0.0:
#                beta_Left[i] = np.max([betaMin,betaMax2])
#                beta_Right[i] = betaMax
#            else:
#                beta_Left[i] = betaMin2
#                beta_Right[i] = np.min([betaMin,betaMax2])
#        
#        self.beta_all = np.concatenate([beta_Left,np.flipud(beta_Right)])
#        self.alpha_all = np.concatenate([alpha_list,np.flipud(alpha_list)])
#        self.psi_bmin = psi_bmin
#        self.psi_bmax = psi_bmax
#
#
#    #               Function to compute the graph alpha vs beta              
#    # =========================================================================
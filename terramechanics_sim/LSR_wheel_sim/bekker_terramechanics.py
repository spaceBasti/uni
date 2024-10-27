# -*- coding: utf-8 -*-
"""
Created on Sat Jan  6 15:25:44 2024

@author: Sebastian Bartel
All formulations in this code are adapted from "NASA White Paper - Terramechanics for LTV Modeling and Simulation" by Zu Qun Li and Lee K. Bingham

Soil Parameters:

n Exponent of np.sinkage 
kc Cohesive modulus 
kφ Frictional modulus 
φ Angle of inertial friction 
c Cohesive strength of soil 
γ Soil weight density 
K Coefficient of soil 
Nq Terzaghi’s bearing capacity factor 
Nc Cohesive bearing capacity factor 
Nγ Density bearing capacity factor 
Kc Cohesive modulus of soil deformation 
Kγ Density modulus of soil deformation 

Wheel Parameters:
    
b Wheel width 
D Wheel diameter 
cf Coefficient of rolling friction 
r Wheel radius 
N Number of wheel grouser 
hg Grouser height

Simulation Input:
    W_wheel Normal force on the wheel 
"""

import numpy as np

def wheel_forces(a, L, h, slope, D, a_x, g, m, mu):
    
    slope = np.deg2rad(slope)
    
    N1 = m*g*((1-a/L)*np.cos(slope)-(h/L+(D/2)/L)*np.sin(slope)-a_x/g)
    N3 = m*g*(a/L*np.cos(slope)+(h/L+(D/2)/L)*np.sin(slope)+a_x/g)  
    T1 = N1/(N1+N3)*(m*g*np.sin(slope)+m*a_x) 
    T3= N3/(N1+N3)*(m*g*np.sin(slope)+m*a_x)      

    #Middle wheels takes up exactly half the load of front and back wheels (only true for identical wheel distance)
    N2 = N1/2 + N3/2
    T2 = T1/2 + T3/2

    #calculate new wheel loading for 6 wheels
    N_tot = N1 + N2 + N3
    N1 = N1 * m*g*np.cos(slope)/N_tot
    N2 = N2 * m*g*np.cos(slope)/N_tot
    N3 = N3 * m*g*np.cos(slope)/N_tot 
    T_tot = T1 + T2 + T3
    T1 = T1 * m*g*np.sin(slope)/T_tot
    T2 = T2 * m*g*np.sin(slope)/T_tot
    T3 = T3 * m*g*np.sin(slope)/T_tot   
    
    return N1, N2, N3, T1, T2, T3

def compression_resistance(n, W_wheel, k_c, b, k_phi, D):
    
    z = (3/(3-n)*W_wheel/((k_c + b*k_phi)*np.sqrt(D)))**(2/(2*n+1))                             #Compression depth [m]
    R_c = (k_c + k_phi*b)*z**(n+1)/(n+1)                                                         
    
    return z, R_c
    
def soil_bearing_capacity(A, c, N_c, gamma, z, N_q, b, N_gamma):
    
    W_s = A*(c*N_c + gamma*z*N_q + 1/2 * gamma * b * N_gamma) 
    
    #Safe Soil Pressure
    p_s_max = W_s/A
    
                                  
    
    return p_s_max
    
def rolling_resistance(W_v, c_f):
    
    R_r = W_v * c_f
    
   
    
    return R_r

def gravitational_resistance(W_v, slope):
    
    R_g = W_v * np.sin(np.deg2rad(slope))
    
    
    
    return R_g

def bulldozing_resistance(b, phi, K_c, K_gamma, z, c, gamma, D):
    
    alpha = np.arccos(1-2*z/D)    
    l_0 = z*np.tan(np.pi/4-phi/2)**2

    R_b = b*np.sin(alpha+phi)/(2*np.sin(alpha)*np.cos(alpha))*(2*z*c*K_c + gamma * z**2*K_gamma) + l_0**3*gamma/3 * (np.pi/2 - phi) + c*l_0**2 * (1 + np.tan(np.pi/4 + phi/2))
    
    #alpha = np.rad2deg(alpha)
    #phi = np.rad2deg(phi)    
    
    #R_b = ((b*np.sin(alpha+phi))/(2*np.sin(alpha)*np.cos(phi)))*(2*c*K_c*z+gamma*K_gamma*z**2) + np.pi*gamma*l_0**2*(90-phi)/540+np.pi*c*l_0**2/180+c*l_0**2*np.tan(45+phi/2)
    
    
    return R_b

def slip(slope):
    
    s = 1/(1+67*np.exp(-0.32*np.deg2rad(slope)))
    
    return s

def max_tractive_force(A, D, b, c, h_g, N , W, phi, K, s, z, l):
    
    N_g = N * l/(D*np.pi)
    H_max = (A*c*(1+2*h_g/b)*N_g+W*np.tan(phi)*(1+0.64*h_g/b*np.arctan(b/h_g))) * (1-K/(s*l) * (1-np.exp(-s*l/K)))
    
    
    #RTS Loco formula
    #H_max =A*c*(1+2*h_g/b)*N_g + W*(1+ (0.64*h_g/b) * np.arctan(b/h_g)) * (1-np.exp(-s * l / K))
    
    
    T_max = H_max * D/2
    
    
    return H_max, T_max


    
    





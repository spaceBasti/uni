# -*- coding: utf-8 -*-
"""
Created on Sat Jan  6 16:31:51 2024

@author: Sebastian Bartel
This code simulates the LSR Locomotion System in the scope of the SEP 23/24 @ FH Aachen

Soil Parameters:

n Exponent of sinkage 
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
mu Friction coefficient

Wheel and Rover Parameters:
    
b Wheel width 
D Wheel diameter 
cf Coefficient of rolling friction 
r Wheel radius 
N_g Number of wheel grouser 
hg Grouser height
rover_mass Mass of rover
W_v Weight of rover
n_wheel Number of wheels
a_x Acceleration in x-direction

Motor and Gear Parameters

i Gear Ratio
eta Efficieny of Gear and Motor
v_nom Nominal translatoric vehicle velocity

"""

import numpy as np
import bekker_terramechanics as bk

print("-------------------------------------------------")
print("Start of Calculation")


#Defining Soil Parameters (Lunar Soil)
n = 1                                                                           #-
k_c = 1400                                                                      #N/m^2
k_phi = 820000                                                                  #N/m^3
phi = 30 * np.pi/180                                                            #rad
c = 170                                                                         #N/m^2
gamma = 2470                                                                    #N/m^3
K = 0.018                                                                       #m
N_q = 32.23                                                                     #-
N_c = 48.09                                                                     #-
N_gamma = 33.27                                                                 #-
K_c = 33.37                                                                     #-
k_gamma = 72.77                                                                 #-
slope = 5                                                                       #°
g = 1.62
mu = 0.2

#Defining Wheel and Rover Parameters
b = 50 * 10**-3                                                                 #m
D = 75 * 10**-3                                                                 #m
c_f = 0.05                                                                      #-
N_g = 12                                                                           #-                                                                           
h_g = 7.5 * 10**-3                                                              #m
rover_mass = 40                                                                 #N
W_v = rover_mass * g                                                            #N
n_wheel = 6                                                                     #-
a = 233.5 * 10**-3                                                              #m
h = 147.5 * 10**-3                                                              #m
l_wheels = 400 * 10**-3                                                         #m
#s = 0.02                                                                        #-

#Defining Gear and Motor Parameters
i = 111                                                                         #-
eta = 0.7 * 0.71                                                                #-
v_nom = 0.2/3.6                                                                 #m/s
t_a = 5                                                                         #s
a_x = v_nom/t_a                                                                 #m/s^2

#Safety Factors
SF_R_c = 1.2                                                                    #Spring
SF_R_b = 1.2                                                                    #Spring
SF_R_g = 1                                                                      #None
SF_R_r = 3                                                                      #Friction


#0. Calculate Wheel Loading dependent on COG position, slope and vehicle acceleration

N1, N2, N3, T1, T2, T3 = bk.wheel_forces(a, l_wheels, h, slope, D, a_x, g, rover_mass, mu) 


#normal forces per wheel
normal_forces = [N1/2, N2/2, N3/2]

#Slip model
s = 0.0353 + 1/(1+67*np.exp(-0.32*slope))
    
    
#Limit Acceleration on Flat Ground
a_limit = g*(1 - a/l_wheels)

#1. Calculate Rolling Resistance 

R_r = bk.rolling_resistance(W_v, c_f)

#2. Calculate Gravitational Resistance

R_g = bk.gravitational_resistance(W_v, slope)



j = 0

z = np.zeros(3)
R_c = np.zeros(3)
l = np.zeros(3)
A = np.zeros(3)
p_s_max = np.zeros(3)
p_s_LSR = np.zeros(3)
H_max = np.zeros(3)
Torque = np.zeros(3)

for N in normal_forces:
    
    #3. Calculate Compression Resistance and Sinkage depth for each wheel

    z[j], R_c[j] = bk.compression_resistance(n, N, k_c, b, k_phi, D)
    
    #4. Calculate contact Area with ground 

    #contact length
    l[j] = D/2 * np.arccos(1-2*z[j]/D)

    #contact Area
    A[j] = l[j]*b

    #5. Calculate maximum Soil Bearing capacity

    p_s_max[j] = bk.soil_bearing_capacity(A[j], c, N_c, gamma, z[j], N_q, b, N_gamma)

    #Current Soil Pressure
    p_s_LSR[j] = N/A[j]
    
    #print("Soil Pressure Wheel %.0f [Pa]: %.3f (Max. Pressure: %.3f" %(j+1) %p_s_LSR[j] %p_s_max[j])
    
    #6. Calculate Tractive Force (Per wheel with grousers)

    #calculate maximum tractive force for each wheel
    #H_max[j], Torque[j] = bk.max_tractive_force(A[j], D, b, c, h_g, N_g, N, phi, K, s, z[j], l[j])
    
    #7. Calculate Bulldozing Resistance (per leading wheel)    
    if j == 0:
        
        R_b = bk.bulldozing_resistance(b, phi, K_c, k_gamma, z[j], c, gamma, D)
    
    j = j + 1


#8. Motor Design

#8.1 DP and Torque per wheel

DP = R_c*SF_R_c + R_b*SF_R_b + R_g/n_wheel + R_r/n_wheel*SF_R_r 
SF = (R_c*SF_R_c + R_b*SF_R_b + R_g/n_wheel + R_r/n_wheel*SF_R_r)/(R_c + R_b + R_g/n_wheel + R_r/n_wheel)
Torque = DP * D/2

#Power per wheel
omega = v_nom/(1-s) * 1/(D/2)
Power = Torque * omega/eta

#Total power
Power_tot = sum(Power) * 2


print("-------------------------------------------------")
print("Torque @ Wheel %.0f (Max. value) @ %.0f °: %.3f" % ((np.argmax(Torque)+1),slope,max(Torque)))
print("Total Power @ %.0f ° and 0.2km/h: %.3f with average SF = %.3f" % (slope, Power_tot, np.average(SF)))



wheel_rpm = omega/(2*np.pi)
motor_rpm = wheel_rpm*i

#Print Values to Command Window

print("-------------------------------------------------")
print("End of Calculation")
print("-------------------------------------------------")































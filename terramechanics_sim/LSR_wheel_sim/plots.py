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
a Longitudinal distance from COG to front axis
h Height of COG measured from center of wheel
l_wheels distance between front and rear wheels

Motor and Gear Parameters

i Gear Ratio
eta Efficieny of Gear and Motor
v_nom Nominal translatoric vehicle velocity

"""

import numpy as np
import bekker_terramechanics as bk
import matplotlib.pyplot as plt


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
slope_range = np.linspace(0,22,1000)                                                                       #°
g = 1.62
mu = 0.2
#s = 0.05

#Defining Wheel and Rover Parameters
b = 50 * 10**-3                                                                 #m
D = 75 * 10**-3                                                                 #m
c_f = 0.05                                                                      #-
N_g = 12                                                                           #-                                                                           
h_g = 7.5 * 10**-3                                                              #m
rover_mass = 39.85                                                                 #kg
W_v = rover_mass * g                                                            #N
n_wheel = 6                                                                     #-
a = 202.1 * 10**-3 #200*10**-3                                                        #m
h = 144.2 * 10**-3                                                              #m
l_wheels = 400 * 10**-3                                                         #m
                                                                        #-

#Defining Gear and Motor Parameters
i = 111                                                                         #-
eta = 0.7 * 0.71                                                                #-
v_nom = 0.2/3.6                                                                 #m/s
t_a = 5                                                                         #s
a_x = v_nom/t_a                                                                 #m/s^2
SF = 2

k = 0
counter = 0

R_g_tot = np.zeros(len(slope_range))
R_c_tot = np.zeros( len(slope_range))
R_b_tot = np.zeros(len(slope_range))
R_r_tot = np.zeros(len(slope_range))
DPPW = np.zeros([3, len(slope_range)])
TPW = np.zeros([3, len(slope_range)])
PPW = np.zeros([3, len(slope_range)])
Power_tot = np.zeros(len(slope_range))
NFPW = np.zeros([3, len(slope_range)])
z = np.zeros([3, len(slope_range)])
DP = np.zeros([3, len(slope_range)])
DP_tot = np.zeros(len(slope_range))
wheel_rpm = np.zeros(len(slope_range))
motor_rpm = np.zeros(len(slope_range))
a_limit = np.zeros(len(slope_range))
t_limit = np.zeros(len(slope_range))


for slope in slope_range:
    
    
    
        
   
    
    #0. Calculate Wheel Loading dependent on COG position, slope and vehicle acceleration

    N1, N2, N3, T1, T2, T3 = bk.wheel_forces(a, l_wheels, h, slope, D, a_x, g, rover_mass, mu)   

    #normal forces per wheel
    normal_forces = [N1/2, N2/2, N3/2]
    shear_forces = [T1/2, T2/2, T3/2]
    
    #Slip model
    s = 0.0353 + 1/(1+67*np.exp(-0.32*slope))
    
    
        
    #Positive Limit Acceleration on depending on slope
    a_limit[k] = g*((1 - a/l_wheels)*np.cos(np.deg2rad(slope))-np.sin(np.deg2rad(slope)))
    t_limit[k] = v_nom/a_limit[k]

    #1. Calculate Gravitational Resistance (for whole vehicle)

    R_g = bk.gravitational_resistance(W_v, slope)
    
    
    
    #Split up Resistance over wheels
    R_g = np.array([R_g/n_wheel, R_g/n_wheel, R_g/n_wheel])



    j = 0

    
    R_c = np.zeros(3)
    R_r = np.zeros(3)
    l = np.zeros(3)
    A = np.zeros(3)
    p_s_max = np.zeros(3)
    p_s_LSR = np.zeros(3)
    Torque = np.zeros(3)
    H_max = np.zeros(3)
    T_max = np.zeros(3)
    
    
    for N in normal_forces:
        
        #2. Calculate Rolling Resistance (for each wheel)

        R_r[j] = bk.rolling_resistance(N, c_f)
        
        #3. Calculate Compression Resistance and Sinkage depth for each wheel

        z[j,k], R_c[j] = bk.compression_resistance(n, N, k_c, b, k_phi, D)
        
        #4. Calculate contact Area with ground 

        #contact length
        l[j] = D/2 * np.arccos(1-2*z[j,k]/D)

        #contact Area
        A[j] = l[j]*b       

        #5. Calculate maximum Soil Bearing capacity

        p_s_max[j] = bk.soil_bearing_capacity(A[j], c, N_c, gamma, z[j,k], N_q, b, N_gamma)

        #Current Soil Pressure
        p_s_LSR[j] = N/A[j]
        
        print(p_s_LSR)

        
        #print("Soil Pressure Wheel %.0f [Pa]: %.3f (Max. Pressure: %.3f" %(j+1) %p_s_LSR[j] %p_s_max[j])
        
        #6. Calculate Tractive Force (Per wheel with grousers) NOT WORKING !!

        #calculate maximum tractive force for each wheel
        #H_max[j], T_max[j] = bk.max_tractive_force(A[j], D, b, c, h_g, N_g, N, phi, K, s, z[j], l[j])
        
        #7. Calculate Bulldozing Resistance (per leading wheel)    
        if j == 0:
            
            R_b = bk.bulldozing_resistance(b, phi, K_c, k_gamma, z[j,k], c, gamma, D)            
            #Bulldozing is zero for Wheels 2-3
            R_b = np.array([R_b, 0, 0])
            
        
        
        
        
        j = j + 1

    #8.1 DP and Torque per wheel NOCHMAL ÜBERPRÜFEN!!!!!
    
    #Drawbar Pull of each wheel
    DP[:,k] = R_c + R_b + R_g/n_wheel + R_r/n_wheel #<- wrong 
    #DP[:,k] = R_c + R_b + R_g + R_r #<- right
    Torque = DP[:,k] * D/2    
    DP_tot[k] = sum(DP[:, k])*2

    #Power per wheel
    omega = (s+1)*v_nom/(D/2)
    #print(omega)
    Power = Torque * omega/eta
    wheel_rpm[k] = omega/(2*np.pi)
    motor_rpm[k] = wheel_rpm[k]*i

    #Total power
    Power_tot[k] = sum(Power) * 2
    
    #Total Resistances
    R_g_tot[k] = sum(R_g)*2
    R_r_tot[k] = sum(R_r)*2
    R_c_tot[k] = sum(R_c)*2
    R_b_tot[k] = R_b[0] * 2    
    
    u=0
    for u in [0, 1, 2]:
        PPW[u, k] = Power[u]
        TPW[u, k] = Torque[u]    
        NFPW[u, k] = normal_forces[u]     
        
    
    
    k = k + 1

#Resistance
fig1 = plt.figure()
plt.grid(True)
plt.plot(slope_range, R_g_tot, label = "Gravitational", color = "black", linestyle = "dashed")
plt.plot(slope_range, R_c_tot, label = "Compaction", color = "black", linestyle = "dotted")
plt.plot(slope_range, R_b_tot, label = "Bulldozing", color = "black", linestyle = "dashdot")
plt.plot(slope_range, R_r_tot, label = "Friction", linestyle = "solid", color = "black")
plt.plot(slope_range, R_g_tot + R_c_tot + R_b_tot + R_r_tot, label = "Total Resistance", color = "grey")
plt.xlim([slope_range[0], 22])
plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22])
plt.ylim([0, 70])
plt.legend(loc='upper left')
plt.xlabel("Slope [°]")
plt.ylabel("Resistance [N]")
plt.savefig('plots/resistance.pdf', format='pdf')

# #Power per Wheel
# fig2 = plt.figure()
# #Wheel 1 (Front Wheel)
# plt.plot(slope_range, PPW[0,:], label = "Front Wheel")
# #Wheel 2 (Middle Wheel)
# plt.plot(slope_range, PPW[1,:], label = "Middle Wheel")
# #Wheel 3 (Rear Wheel)
# plt.plot(slope_range, PPW[2,:], label = "Rear Wheel")
# plt.grid(True)
# plt.xlabel("Slope [°]")
# plt.ylabel("Power [W]")
# plt.legend()
# plt.savefig("plots/power_per_wheel.pdf", format='pdf')

#Total Power with SF
ind_mean_slope = np.where(np.isclose(slope_range, 5, atol=0.01))
ind_max_slope = np.where(np.isclose(slope_range, 11, atol=0.01))
ind_ingress_slope = np.where(np.isclose(slope_range, 20, atol=0.01))

fig3 = plt.figure()
plt.plot(slope_range, SF*Power_tot, linestyle = "-", color = "black")
plt.annotate("P=%.2f W" %(SF*Power_tot[ind_mean_slope]), [slope_range[ind_mean_slope]+0.5, SF*Power_tot[ind_mean_slope]-0.5])
plt.annotate("P=%.2f W" %(SF*Power_tot[ind_max_slope]), [slope_range[ind_max_slope]+0.5, SF*Power_tot[ind_max_slope]-0.5])
print(SF*Power_tot[ind_mean_slope])
print(SF*Power_tot[ind_max_slope])
print(SF*Power_tot[ind_ingress_slope])
plt.vlines(5,6, 9)
plt.vlines(11,6, 11)
plt.grid(True)
plt.xlabel("Slope [°]")
plt.ylabel("Total Power (incl. SF = %.2f) [W]" %SF)
#plt.xlim([slope_range[0], 22])
plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22])
#plt.ylim([6, 15])
plt.savefig('plots/total_power.pdf', format='pdf')

#Wheel Loads
fig4 = plt.figure()
#Wheel 1 (Front Wheel)
plt.plot(slope_range, NFPW[0,:], label = "Front Wheel")
#Wheel 2 (Middle Wheel)
plt.plot(slope_range, NFPW[1,:], label = "Middle Wheel")
#Wheel 3 (Rear Wheel)
plt.plot(slope_range, NFPW[2,:], label = "Rear Wheel")
plt.grid(True)
plt.xlabel("Slope [°]")
plt.ylabel("Normal Load [N]")
plt.legend()
plt.title("Normal Loads")
plt.savefig('plots/wheel_loads.pdf', format='pdf')

#Sinkage Depth
fig5 = plt.figure()
#Wheel 1 (Front Wheel)
plt.plot(slope_range, z[0,:]*1000, label = "Front Wheel", linestyle = "solid", color = "black")
#Wheel 2 (Middle Wheel)
plt.plot(slope_range, z[1,:]*1000, label = "Middle Wheel", linestyle = "dashed", color = "black")
#Wheel 3 (Rear Wheel)
plt.plot(slope_range, z[2,:]*1000, label = "Rear Wheel", linestyle = "dashdot", color = "black")
plt.grid(True)
plt.xlim([slope_range[0], 22])
plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22])
plt.ylim([6, 18])
plt.xlabel("Slope [°]")
plt.ylabel("Sinkage Depth [mm]")
plt.legend()

plt.savefig('plots/sinkage_depth.pdf', format='pdf')

#Drawbar Pull for each wheel
fig6 = plt.figure()
#Wheel 1 (Front Wheel)
plt.plot(slope_range, DP[0,:], label = "Front Wheel", linestyle = "solid", color = "black")
#Wheel 2 (Middle Wheel)
plt.plot(slope_range, DP[1,:], label = "Middle Wheel", linestyle = "dashed", color = "black")
#Wheel 3 (Rear Wheel)
plt.plot(slope_range, DP[2,:], label = "Rear Wheel", linestyle = "dashdot", color = "black")
plt.grid(True)
plt.xlim([slope_range[0], 22])
plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22])
#plt.ylim([6, 18])
plt.xlabel("Slope [°]")
plt.ylabel("DP per wheel [N]")
plt.legend()

plt.savefig('plots/drawbar_pull.pdf', format='pdf')

#Total Drawbar Pull 
fig7 = plt.figure()
#Wheel 1 (Front Wheel)
plt.plot(slope_range, DP_tot, linestyle = "solid", color = "black")
plt.grid(True)
plt.xlim([slope_range[0], 22])
plt.xticks([0,2,4,6,8,10,12,14,16,18,20,22])
#plt.ylim([6, 18])
plt.xlabel("Slope [°]")
plt.ylabel("Total DP [N]")


plt.savefig('plots/drawbar_pull_total.pdf', format='pdf')

#RPM 
fig8, ax8 = plt.subplots()
ax8.plot(TPW[0,:]/i*1000, motor_rpm*60, linestyle = "solid", color = "black", label = "Front Wheel", )
ax8.plot(TPW[1,:]/i*1000, motor_rpm*60, linestyle = "dashed", color = "black", label = "Middle Wheel")
ax8.plot(TPW[2,:]/i*1000, motor_rpm*60, linestyle = "dashdot", color = "black", label = "Rear Wheel")
plt.grid(True)
plt.xlim([1.25, 4])
plt.xticks([1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5, 3.75, 4, 4.25])
plt.ylim([1600, 3200])
ax8.set_xlabel("Front motor torque [mNm]")
ax8.set_ylabel("Motor speed [RPM]")
plt.legend()

plt.savefig('plots/rpm_torque.pdf', format='pdf')



# -*- coding: utf-8 -*-

===============================================================================
Copyright 2025 Akshay Harlalka 

Thermal Model for Freeze Preventive Vaccine Carriers
"""


import numpy as np
import matplotlib.pyplot as plt

rho1 = 200.0 # density of polyurethane in kg/m3 
rho2 = 920.0 # density of ice in kg / m^3
rho3 = 1000 # density of 'buffer layer 1' in kg / m^3
rho4 = 200 # density of 'buffer layer 2' in kg / m^3
rho5 = 1000 # density of vaccine in kg / m^3

k1 = 0.026 # thermal conductivity of polyurethane in W / (m*K)
k2 = 2.22 # thermal conductivity of ice in W / (m*K)
k3 = 0.6 # thermal conductivity of buffer layer 1 in W / (m*K)
k4 = 0.026 # thermal conductivity of buffer layer 2 in W / (m*K)
k5 = 0.6 # thermal conductivity of vaccine in W / (m*K)

c1 = 1400.0 # specific heat capacity of polyurethane in J / (kg*K)
c2 = 2040.0 # specific heat capacity of ice in J / (kg*K)
c3 = 4182 # specific heat capacity of thermal buffer layer 1 in J / (kg*K)
c4= 1400.0 # specific heat capacity of thermal buffer layer 2 in J / (kg*K)
c5 = 4182.0 # specific heat capacity of vaccine in J / (kg*K)
h1 = 25.0 # convective heat transfer coefficient in W / (m^2 * K)

T_initial1 = 316.0 # initial temperature of PU in Kelvin
T_initial2 = 248.0 # initial temperature of ice in Kelvin
T_initial3 = 283 #initial temperature of the thermal buffer layer 1 in kelvin
T_initial4 = 283 #initial temperature of the thermal layer 2 in kelvin
T_initial5 = 279 #initial temperature of the vaccine in kelvin
T_inf1 = 316.0 # ambient temperature in Kelvin

#Thickness of different control volumes
dx1=0.0425 
dx2=0.033
dx3=0.012 
dx4=0.012 
dx5=0.016 

total_time = 3600*30 # total duration of simulation in seconds
nsteps = 3600 # number of timesteps
dt = total_time/nsteps # duration of timestep in seconds

# initialize volume element coordinates and time samples
timesamps = np.linspace(0, dt*nsteps, nsteps+1)

# initialize a big 2D array to store temperature values

T=np.zeros([5,nsteps+1])
T[0,0]=T_initial1
T[1,0]=T_initial2
T[2,0]=T_initial3
T[3,0]=T_initial4
T[4,0]=T_initial5


#Heat Balance Equations for each control volume
for j in range(len(timesamps)-1):
    T[0,j+1]= T[0,j]+ ((h1*dt*(T_inf1-T[0,j]))/(rho1*c1*dx1)) + (k1*dt*(T[1,j]-T[0,j])/(rho1*c1*(dx1*dx1)))
    T[1,j+1]=T[1,j]+ ((k1*(T[0,j]-T[1,j])*dt)/(dx1*rho2*c2*dx2)) + ((k2*(T[2,j]-T[1,j])*dt)/(dx2*rho2*c2*dx2))
    T[2,j+1]=T[2,j]+ ((k2*(T[1,j]-T[2,j])*dt)/(dx2*rho3*c3*dx3)) + ((k3*(T[3,j]-T[2,j])*dt)/(dx3*rho3*c3*dx3))
    T[3,j+1]=T[3,j]+ ((k3*(T[2,j]-T[3,j])*dt)/(dx3*rho4*c4*dx4)) + ((k4*(T[4,j]-T[3,j])*dt)/(dx4*rho4*c4*dx4))
    T[4,j+1]=T[4,j] + (k4*dt*(T[3,j]-T[4,j])/(rho5*c5*dx4*dx5))

# this plots the temperature vs time plots
plt.figure(1)
T1=T[0,:]-273
T2=T[1,:]-273
T3=T[2,:]-273
T4=T[3,:]-273
T5=T[4,:]-273
plt.plot(np.linspace(0,nsteps,num=nsteps+1)*dt,T5)
plt.xlabel('Time (secs)')
plt.ylabel('Temperature (C)')

plt.show()




# end of file

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp
from vpython import *
from mpl_toolkits.mplot3d import Axes3D

# Parameters
global Na
global Nb
global N
Na = 800
Nb = 200
N = Na + Nb
sigmaAA, sigmaAB, sigmaBB = 1.0, 0.8, 0.88
sigmaAA6, sigmaAB6, sigmaBB6= sigmaAA**6,sigmaAB**6, sigmaBB**6
sigmaAA12, sigmaAB12, sigmaBB12 = sigmaAA6 **2, sigmaAB6**2, sigmaBB6**2
epsAA, epsAB, epsBB = 1.0, 1.5, 0.5
m = 1 # mA = mB
t0 = 0.0
dt = 0.05
dt_sq = dt**2
temp = 1.0

# potential cutoff
rcutAA = 2.5*sigmaAA
rcutAB = 2.5*sigmaAB
rcutBB = 2.5*sigmaBB

rcutAA_sq = rcutAA**2
rcutAB_sq = rcutAB**2
rcutBB_sq = rcutBB**2

x = np.zeros(N, float)
y= np.zeros(N, float)
z= np.zeros(N, float)

x,y,z = np.loadtxt('initposvel', dtype = float ,unpack = True)

# Randomly distribute the particles in the box
np.random.seed(15)
L = 9.4 # Linear Dimension of the box
Ldiv2 = L/2.0
x,y,z = np.random.uniform(low=0.0, high=L, size=(3,N))

# Plot the initial configuration
#3d plot
fig3d = plt.figure()
fax = fig3d.add_subplot(111, projection='3d')
fax.scatter(x[:Na], y[:Na], z[:Na], c='r', marker='o', label='A')
fax.scatter(x[Na:], y[Na:], z[Na:], c='b', marker='o', label='B')
fax.set_xlabel('X')
fax.set_ylabel('Y')
fax.set_zlabel('Z')
plt.title('Initial Configuration')
plt.legend()
plt.show()

# Velocities from Maxwell-Boltzmann Distribution
def maxwellboltzmanvel(temp):
    global vx
    global vy
    global vz
    nopart=len(vx)
    sigma = np.sqrt(temp) # sqrt. (kT/m) = sqrt. (T)
    vx = np.random.normal(0.0, sigma, nopart)
    vy = np.random.normal(0.0, sigma, nopart)
    vz = np.random.normal(0.0, sigma, nopart)
    
    # Checking that center of mass doesn't drift
    vx -= sum(vx)/float(nopart)
    vy -= sum(vy)/float(nopart)
    vz -= sum(vz)/float(nopart)
    
    # Checking Temperature
    scalefactor = np.sqrt(3.0*temp*nopart/sum(vx**2 + vy**2 + vz**2))
    vx *= scalefactor
    vy *= scalefactor
    vz *= scalefactor

# arrow(pos=vector(tx,ty,tz), axis=vector(tvx,tvy,tvz), color=color.green)

with open('./trajectory_data.txt', 'w') as traj_file:
        traj_file.write("time, particle, x, y, z, vx, vy, vz\n") # Header for clarity
# Accelerations
def acceleration(x,y,z, Na, N):
    global L
    global Ldiv2
    
    ax = sp.zeros(N)
    ay = sp.zeros(N)
    az = sp.zeros(N)
    
    # AA interactions
    for i in range(N-1):
        xi, yi, zi = x[i], y[i], z[i]
        
        for j in range(i+1, N):
            xij=xi-x[j]
            yij=yi-y[j]
            zij=zi-z[j]
            
            # Minimum Image Convention
            if xij > Ldiv2:xij-=L
            elif xij<-Ldiv2:xij+=L
            
            if yij > Ldiv2:yij-=L
            elif yij<-Ldiv2:yij+=L
            
            if zij > Ldiv2:zij-=L
            elif zij<-Ldiv2:zij+=L
            
            rijto2=xij**2+yij**2+zij**2
            if(rijto2 < rcutto2):
                onedivrijto2 = 1.0/rijto2
                fmagtmp = eps * (sigmato12*onedivrijto2**7 - 0.5*sigmato6*onedivrijto2**4)
                ax[i] += fmagtmp*xij
                ax[j] -= fmagtmp*xij
                ay[i] += fmagtmp*yij
                ay[j] -= fmagtmp*yij
                az[i] += fmagtmp*zij
                az[j] -= fmagtmp*zij
                
                if i < Na:
                    if j < Na: # AA
                        rcutto2 = rcutAA_sq
                        sigmato12 = sigmaAA12
                        sigmato6 = sigmaAA6
                        eps = epsAA
                    else: #AB
                        rcutto2 = rcutAB_sq
                        sigmato12 = sigmaAB12
                        sigmato6 = sigmaAB6
                        eps = epsAB
                else: #BB
                    rcutto2 = rcutBB_sq
                    sigmato12 = sigmaBB12
                    sigmato6 = sigmaBB6
                    eps = epsBB
        
        # AB interactions
        for i in range(Na):
            xi, yi, zi = x[i], y[i], z[i]
            for j in range(Na, N):
                # Minimum image convention
                if xij > Ldiv2:xij-=L
                elif xij<-Ldiv2:xij+=L
                
                if yij > Ldiv2:yij-=L
                elif yij<-Ldiv2:yij+=L
                
                if zij > Ldiv2:zij-=L
                elif zij<-Ldiv2:zij+=L
                
                rijto2 = xij*xij + yij*yij + zij*zij
                if (rijto2 < rcutAB_sq):
                    onedivrijto2 = 1.0/rijto2
                    fmagtmp = epsAB*(sigmaAB12*onedivrijto2**7 - 0.5*sigmaAB6*onedivrijto2**4)
                    ax[i] += fmagtmp*xij
                    ax[j] -= fmagtmp*xij
                    ay[i] += fmagtmp*yij
                    ay[j] -= fmagtmp*yij
                    az[i] += fmagtmp*zij
                    az[j] -= fmagtmp*zij
                    
                    if i < Na:
                        if j < Na: # AA
                            rcutto2 = rcutAA_sq
                            sigmato12 = sigmaAA12
                            sigmato6 = sigmaAA6
                            eps = epsAA
                        else: #AB
                            rcutto2 = rcutAB_sq
                            sigmato12 = sigmaAB12
                            sigmato6 = sigmaAB6
                            eps = epsAB
                    else: #BB
                        rcutto2 = rcutBB_sq
                        sigmato12 = sigmaBB12
                        sigmato6 = sigmaBB6
                        eps = epsBB
        
        # BB interactions
        for i in range(Na, N-1):
            xi, yi, zi = x[i], y[i], z[i]
            
            for j in range(i+1, N):
                xij=xi-x[j]
                yij=yi-y[j]
                zij=zi-z[j]
                
                # Minimum Image Convention
                if xij > Ldiv2:xij-=L
                elif xij<-Ldiv2:xij+=L
                
                if yij > Ldiv2:yij-=L
                elif yij<-Ldiv2:yij+=L
                
                if zij > Ldiv2:zij-=L
                elif zij<-Ldiv2:zij+=L
                
                rijto2=xij**2+yij**2+zij**2
                if(rijto2 < rcutto2):
                    onedivrijto2 = 1.0/rijto2
                    fmagtmp = eps * (sigmato12*onedivrijto2**7 - 0.5*sigmato6*onedivrijto2**4)
                    ax[i] += fmagtmp*xij
                    ax[j] -= fmagtmp*xij
                    ay[i] += fmagtmp*yij
                    ay[j] -= fmagtmp*yij
                    az[i] += fmagtmp*zij
                    az[j] -= fmagtmp*zij
                    
                    if i < Na:
                        if j < Na: # AA
                            rcutto2 = rcutAA_sq
                            sigmato12 = sigmaAA12
                            sigmato6 = sigmaAA6
                            eps = epsAA
                        else: #AB
                            rcutto2 = rcutAB_sq
                            sigmato12 = sigmaAB12
                            sigmato6 = sigmaAB6
                            eps = epsAB
                    else: #BB
                        rcutto2 = rcutBB_sq
                        sigmato12 = sigmaBB12
                        sigmato6 = sigmaBB6
                        eps = epsBB

    return 48*ax, 48*ay, 48*az

print("Trajectory data saved to trajectory_data.txt")
ax, ay, az = acceleration(x, y, z)

# MD Simulation via Velocity Verlet
for tstep in range(1, N+1):
    # update positions
    x += vx*dt + 0.5*ax*dt_sq
    y += vy*dt + 0.5*ay*dt_sq
    z += vz*dt + 0.5*az*dt_sq
    #periodic boundary conditions
    for i in range(N):
        if x[i]>L: x[i]-=L
        elif x[i]<=0: x[i]+= L
        if y[i]>L: y[i]+= L
        elif y[i] <= 0: y[i] -=L
        if z[i] > L: z[i]-=L
        elif z[i] <= L: z[i]+=L
    
    #update velocities
    vx += 0.5*ax*dt
    vy += 0.5*ay*dt 
    vz += 0.5*az*dt
    ax, ay, az = acceleration(x, y, z)
    vx += 0.5*ax*dt
    vy += 0.5*ay*dt 
    vz += 0.5*az*dt
    
    current_time = tstep * dt
    for i in range(N):
        with open('trajectory_data.txt', 'w') as traj_file:
            traj_file.write("time, i, x, y, z, vx, vy, vz, ax, ay, az\n") # Header for clarity
        traj_file.write(f"{current_time},{i},{x[i]},{y[i]},{z[i]},{vx[i]},{vy[i]},{vz[i]}, {ax[i]}, {ay[i]}, {az[i]}\n")
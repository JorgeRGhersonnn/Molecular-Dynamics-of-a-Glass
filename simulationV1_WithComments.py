#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp
# from vpython import *
from mpl_toolkits.mplot3d import Axes3D

global Na
global Nb
global N
Na = 800
Nb = 200
N = Na + Nb

x = np.zeros(N, float)
y= np.zeros(N, float)
z= np.zeros(N, float)

#KVL: please note that initposvel actually also has vx,vy,vz
#         I see that you have the Boltzmann-distributed vel, but this
#         initposvel allow you to first get your NVE simulation working
#         without having to carefully create your own initial config
x,y,z = np.loadtxt('initposvel', dtype = float ,unpack = True)

# Randomly distribute the particles in the box
np.random.seed(15)
L = 9.4 # Linear Dimension of the box
Ldiv2 = L/2.0
#KVL: I highly recommend against that for now, because then some of your particles might end up very close and so you might get huge accelerations. If you want to do MD at finite temperature, then that is a big problem, plus there is no need for it, because you have     initposvel
x,y,z = np.random.uniform(low=0.0, high=L, size=(3,N))

#KVL: okay, sorry, in class when I looked at it, I misunderstood what you are doing here. I still recommend that you start instead with initpos, because the cubic lattice might be difficult to get melted. Unless, you could do simulations starting with e.g. initposvel and also another simulation(s) starting with lattice and you could compare results, that might be kind of cool.   Here just notice that your x,y,z   values will be the ones which you had last in program, so you started with cubic lattice. If you don't want this, you need to comment out the other initializations. Huge blocks you can comment out with '''   then text block and end with    '''
# Cubic Lattice
nsitesx = int(round(pow(N,1.0/3.0)))
dsitesx = L/float(nsitesx)

for ni in range(nsitesx):
    tmpz = (0.5+ni)*dsitesx
    for nj in range(nsitesx):
        tmpy = (0.5 + nj)*dsitesx
        for nk in range(nsitesx):
            tmpx = (0.5 + nk)*dsitesx
            i = nk+nj*nsitesx+ni*nsitesx**2
            x[i] = tmpx
            y[i] = tmpy
            z[i] = tmpz

# Randomly shuffle the particles
for i in range(Na, N):
    j = np.random.randint(Na)
    x[i], x[j] = x[j], x[i]
    y[i], y[j] = y[j], y[i]
    z[i], z[j] = z[j], z[i]

np.savetx('initposcheck', (sp.transpose(sp.vstack(x,y,z))))
    
# Plot the initial configuration
#2d plot
fig = plt.figure()
plt.scatter(x[:Na], z[:Na], s=150, c='r', marker='o', label='A')
plt.scatter(x[Na:], z[Na:], s=70, c='b', marker='o', label='B')
plt.xlim(0,L)
plt.xlabel('$x$')
plt.ylabel('$z$')
plt.title('Initial Configuration')
plt.legend()
plt.show()

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

# Accelerations
def acceleration(x,y,z):
    global L
    global Ldiv2
    
    ax = sp.zeros(N)
    ay = sp.zeros(N)
    az = sp.zeros(N)
    #KVL: in AJP paper bottom page 408 it explains that you might want to use three separate loops over i and j, one for AA, one for AB, and one for BB, 
    for i in range(N-1):
        xi=x[i]
        yi=y[i]
        zi=z[i]
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
            #KVL: continue here

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
# from vpython import *

# Accelerations
def acceleration(x,y,z):
    '''
    Inputs:
    - x-axis position. 
    - y-axis position.
    - z-axis position.
    
    Outputs: 
    - x-axis acceleration 
    - y-axis acceleration
    - z-axis acceleration
    '''
    global L
    global Ldiv2
    global Na,Nb,N
    global epsAA,epsAB,epsBB
    global sigmaAA12,sigmaAB12,sigmaBB12
    global sigmaAA6,sigmaAB6,sigmaBB6
    global rcutAA_sq,rcutAB_sq,rcutBB_sq
    
    # Initialize accelerations
    ax = np.zeros(N)
    ay = np.zeros(N)
    az = np.zeros(N)
    
    # A-A Interactions
    for i in range(Na-1):
        xi = x[i]
        yi = y[i]
        zi = z[i]        
        for j in range(i+1, N):
            xij = xi - x[j]
            yij = yi - y[j]
            zij = zi - z[j]
            # minimum-image
            if xij > Ldiv2: 
                xij -= L
            elif xij < -Ldiv2: 
                xij += L
            if yij > Ldiv2: 
                yij -= L
            elif yij < -Ldiv2: 
                yij += L
            if zij > Ldiv2: 
                zij -= L
            elif zij < -Ldiv2: 
                zij += L    
                
            r2 = xij*xij + yij*yij + zij*zij
            if(r2 < rcutAA_sq):
                inv_r2 = 1.0 / r2
                # Direct force components:
                # F_x = 48 * eps * [ sigma**12 dx/r**14 − 0.5*(sigma**6 dx/r**8)]
                fmagtmp = epsAA*(sigmaAA12*inv_r2**7 - 0.5*sigmaAA6*inv_r2**4)
                ax[i] += fmagtmp*xij
                ax[j] -= fmagtmp*xij
                ay[i] += fmagtmp*yij
                ay[j] -= fmagtmp*yij
                az[i] += fmagtmp*zij
                az[j] -= fmagtmp*zij
                
    # A-B interactions
    for i in range(Na): 
        xi = x[i]
        yi = y[i]
        zi = z[i]        
        for j in range(Na, N):
            xij = xi = x[j]
            yij = yi - y[j]
            zij = zi - z[j]
            # minimum-image
            if xij > Ldiv2: 
                xij -= L
            elif xij < -Ldiv2: 
                xij += L
            if yij > Ldiv2: 
                yij -= L
            elif yij < -Ldiv2: 
                yij += L
            if zij > Ldiv2: 
                zij -= L
            elif zij < -Ldiv2: 
                zij += L    
            
            r2 = xij*xij + yij*yij + zij*zij
            if(r2 < rcutAB_sq):
                inv_r2 = 1.0 / r2
                # Direct force components:
                # F_x = 48 * eps * [ sigma**12 dx/r**14 − 0.5*(sigma**6 dx/r**8)]
                fmagtmp = epsAB*(sigmaAB12*inv_r2**7 - 0.5*sigmaAB6*inv_r2**4)
                ax[i] += fmagtmp*xij
                ax[j] -= fmagtmp*xij
                ay[i] += fmagtmp*yij
                ay[j] -= fmagtmp*yij
                az[i] += fmagtmp*zij
                az[j] -= fmagtmp*zij
            
    # B-B Interactions
    for i in range(Na, N-1): 
        xi = x[i]
        yi = y[i]
        zi = z[i]        
        for j in range(Na, N-1):
            xij = xi - x[j]
            yij = yi - y[j]
            zij = zi - z[j]
            # minimum-image
            if xij > Ldiv2: 
                xij -= L
            elif xij < -Ldiv2: 
                xij += L
            if yij > Ldiv2: 
                yij -= L
            elif yij < -Ldiv2: 
                yij += L
            if zij > Ldiv2: 
                zij -= L
            elif zij < -Ldiv2: 
                zij += L    
            
            r2 = xij*xij + yij*yij + zij*zij
            if(r2 < rcutBB_sq):
                inv_r2 = 1.0 / r2
                # Direct force components:
                # F_x = 48 * eps * [ sigma**12 dx/r**14 − 0.5*(sigma**6 dx/r**8)]
                fmagtmp = epsBB*(sigmaBB12*inv_r2**7 - 0.5*sigmaBB6*inv_r2**4)
                ax[i] += fmagtmp*xij
                ax[j] -= fmagtmp*xij
                ay[i] += fmagtmp*yij
                ay[j] -= fmagtmp*yij
                az[i] += fmagtmp*zij
                az[j] -= fmagtmp*zij
                
    return 48*ax, 48*ay, 48*az
            
# Kob-Andersen potential function 
def potential(x,y,z):
    '''
    Inputs:
    - x-axis position.
    - y-axis position.
    - z-axis position
    
    Output:
    - Potential 
    '''
    global L 
    global Na, Nb, N
    
    Vtot = 0.0
    
    # A-A interactions
    for i in range(Na-1):
        xi=x[i]
        yi=y[i]
        zi=z[i]
        for j in range(i+1,Na):
            xij=xi-x[j]
            yij=yi-y[j]
            zij=zi-z[j]
            #minimum image convention
            if xij > Ldiv2: xij -= L
            elif xij < - Ldiv2: xij  += L
            if yij > Ldiv2: yij -= L
            elif yij < - Ldiv2: yij  += L
            if zij > Ldiv2: zij -= L
            elif zij < - Ldiv2: zij  += L

            rijto2 = xij*xij + yij*yij + zij*zij
            if(rijto2 < rcutAA_sq):
                onedivrijto2 = 1.0/rijto2
                Vtot += epsAA*(sigmaAA12*onedivrijto2**6
                            - sigmaAA6*onedivrijto2**3)-VcutAA 
                
    # AB interactions
    for i in range(Na):
        xi=x[i]
        yi=y[i]
        zi=z[i]
        for j in range(Na,N):
            xij=xi-x[j]
            yij=yi-y[j]
            zij=zi-z[j]
            #minimum image convention
            if xij > Ldiv2: xij -= L
            elif xij < - Ldiv2: xij  += L
            if yij > Ldiv2: yij -= L
            elif yij < - Ldiv2: yij  += L
            if zij > Ldiv2: zij -= L
            elif zij < - Ldiv2: zij  += L

            rijto2 = xij*xij + yij*yij + zij*zij
            if(rijto2 < rcutAB_sq):
                onedivrijto2 = 1.0/rijto2
                Vtot += epsAB*(sigmaAB12*onedivrijto2**6
                            - sigmaAB6*onedivrijto2**3)-VcutAB

    # BB interactions
    for i in range(Na,N-1):
        xi=x[i]
        yi=y[i]
        zi=z[i]
        for j in range(i+1,N):
            xij=xi-x[j]
            yij=yi-y[j]
            zij=zi-z[j]
            #minimum image convention
            if xij > Ldiv2: xij -= L
            elif xij < - Ldiv2: xij  += L
            if yij > Ldiv2: yij -= L
            elif yij < - Ldiv2: yij  += L
            if zij > Ldiv2: zij -= L
            elif zij < - Ldiv2: zij  += L

            rijto2 = xij*xij + yij*yij + zij*zij
            if(rijto2 < rcutBB_sq):
                onedivrijto2 = 1.0/rijto2
                Vtot += epsBB*(sigmaBB12*onedivrijto2**6
                            - sigmaBB6*onedivrijto2**3)-VcutBB

    return 4.0*Vtot/float(N)

# radial distribution measurement
def histmeas(x,y,z):
    global histgofrAA,histgofrBB,histgofrAB
    global L
    global Ldiv2
    global Na, Nb, N
    global grnbinmax
    global grdelta
    for i in range(0,N-1):
        xi=x[i]
        yi=y[i]
        zi=z[i]
        for j in range(i+1,N):
            xij=xi-x[j]
            yij=yi-y[j]
            zij=zi-z[j]
            #minimum image convention
            if xij > Ldiv2: 
                xij -= L
            elif xij < - Ldiv2: 
                xij  += L
            if yij > Ldiv2: 
                yij -= L
            elif yij < - Ldiv2: 
                yij  += L
            if zij > Ldiv2: 
                zij -= L
            elif zij < - Ldiv2: 
                zij  += L

            rijto2 = xij*xij + yij*yij + zij*zij
            
            # Checker
            if (not np.isfinite(rijto2) or rijto2 <= 0):
                continue
            
            rij= np.sqrt(rijto2)
            grbin=int(rij/grdelta)
            if(grbin < grnbinmax):
                if(i < Na):
                    if (j < Na): #AA
                        histgofrAA[grbin] += 2
                    elif(): 
                        histgofrAB[grbin] += 1
                    else: #BB
                        histgofrBB[grbin] += 2

# Maxwell Boltzmann distribution
def maxwellboltzmannvel(temp):
    global vx
    global vy
    global vz
    nopart=len(vx)
    sigma=np.sqrt(temp) #sqrt(kT/m)
    vx=np.random.normal(0.0,sigma,nopart)
    vy=np.random.normal(0.0,sigma,nopart)
    vz=np.random.normal(0.0,sigma,nopart)
    #  make sure that center of mass does not drift
    vx -= sum(vx)/float(nopart)
    vy -= sum(vy)/float(nopart)
    vz -= sum(vz)/float(nopart)
    #  make sure that temperature is exactly wanted temperature
    scalefactor = np.sqrt(3.0*temp*nopart/sum(vx*vx+vy*vy+vz*vz))
    vx *= scalefactor
    vy *= scalefactor
    vz *= scalefactor

# ——————————— Initialize ———————————
global Na
global Nb
global N
# Parameters
Na      = 800               # # of A particles
Nb      = 200               # # of B particles
N       = Na + Nb           # total
L       = 9.4               # box length
Ldiv2   = L/2.0
dt      = 0.005             # time step
dtto2   = dt**2
nstepMD = 1000              # # of MD steps
rcutfactor = 2.5

# KA LJ parameters (units: epsAA=1, sigmaAA=1)
sigmaAA=1.0
sigmaAB=0.8
sigmaBB=0.88
epsAA=1.0
epsAB=1.5
epsBB=0.5

sigmaAA12 = sigmaAA**12
sigmaAB12 = sigmaAB**12
sigmaBB12 = sigmaBB**12
sigmaAA6 = sigmaAA**6
sigmaAB6 = sigmaAB**6
sigmaBB6 = sigmaBB**6

# radius cutoffs
rcutAA_sq = (rcutfactor*sigmaAA)**2
rcutAB_sq = (rcutfactor*sigmaAB)**2
rcutBB_sq = (rcutfactor*sigmaBB)**2

# Potential Cutoffs
VcutAA = epsAA*(sigmaAA12/(rcutAA_sq**6)-sigmaAA6/(rcutAA_sq**3))
VcutAB = epsAB*(sigmaAB12/(rcutAB_sq**6)-sigmaAB6/(rcutAB_sq**3))
VcutBB = epsBB*(sigmaBB12/(rcutBB_sq**6)-sigmaBB6/(rcutBB_sq**3))

# for g(r) histogram
grdelta = 0.1
grnbinmax = int(Ldiv2/grdelta)
histgofrAA = np.zeros(grnbinmax,int)
histgofrAB = np.zeros(grnbinmax,int)
histgofrBB = np.zeros(grnbinmax,int)
nstepgofr = 25

np.random.seed(15)

x = np.zeros(N,float)
y = np.zeros(N,float)
z = np.zeros(N,float)
vx = np.zeros(N,float)
vy = np.zeros(N,float)
vz = np.zeros(N,float)

# load initposvel: x,y,z, vx,vy,vz
x, y, z, vx, vy, vz = np.loadtxt('initposvel', dtype='float', unpack=True)

# Boltzmann Dist. parameters
nstepBoltz = 1
temp = 0.5

# * 2d-scatter plot (for data points)
plt.figure()
plt.scatter(x[:Na],y[:Na],s=150,color='blue')
plt.scatter(x[Na:],y[Na:],s=70,color='red')
plt.xlabel('$v_x$')
plt.ylabel('$v_z$')
plt.show()

# * 3d-Vpython figure
'''
for i in range(N):
  tx=x[i]
  ty=y[i]
  tz=z[i]
  tvx = vx[i]
  tvy = vy[i]
  tvz = vz[i]
  if i < Na :
    sphere(pos=vector(tx,ty,tz),radius=0.5,color=color.blue)
  else:
    sphere(pos=vector(tx,ty,tz),radius=0.2,color=color.red)
  arrow(pos=vector(tx,ty,tz),axis=vector(tvx,tvy,tvz),color=color.green)
# With middle mouse button zoom in & out; with right mouse button rotate
'''
# initialize accelerations
ax, ay, az = acceleration(x,y,z)

# check accelerations
np.savetxt('initacccheck', (np.transpose(np.vstack((ax,ay,az)))))

temp_array = np.zeros(nstepMD,float)

# ————————— Velocity‐Verlet MD —————————
for tstep in range(1, nstepMD+1):
    # update positions
    x += vx*dt + 0.5*ax*dtto2
    y += vy*dt + 0.5*ay*dtto2
    z += vz*dt + 0.5*az*dtto2
    
    # periodic boundary conditions:
    for i in range(N):
        if x[i] > L: 
            x[i] -= L
        elif x[i] <= 0:  
            x[i] += L
        if y[i] > L: 
            y[i] -= L
        elif y[i] <= 0: 
            y[i] += L
        if z[i] > L: 
            z[i] -= L
        elif z[i] <= 0:  
            z[i] += L

    # update velocities
    vx += 0.5*ax*dt
    vy += 0.5*ay*dt
    vz += 0.5*az*dt
    ax,ay,az = acceleration(x,y,z)
    vx += 0.5*ax*dt
    vy += 0.5*ay*dt
    vz += 0.5*az*dt

    # Temperature bath
    if(tstep % nstepBoltz) == 0:
        maxwellboltzmannvel(temp)

    # Determine KE and Temperature
    KE = 0.5*sum(vx*vx + vy*vy + vz*vz)/float(N)
    temp_meas = KE*(2.0/3.0)
    temp_array[tstep - 1] = temp_meas  
    
    # Measure RDF
    if(tstep % nstepgofr) == 0:
        histmeas(x,y,z)
    
# Print r and RDF for AA, AB, BB
fileoutgofr = open("gofrAABBAB.data", "w")
meas = int(nstepMD/nstepgofr)
for grbin in range(grnbinmax):
    rinner = grbin*grdelta
    router = rinner + grdelta
    shellvol = (4.0*np.pi/3.0)*(router**3 - rinner**3)
    gofrAA = (L**3/(Na*(Na-1)))*histgofrAA[grbin]/(shellvol*meas)
    gofrBB = (L**3/(Nb*(Nb-1)))*histgofrBB[grbin]/(shellvol*meas)
    gofrAB = (L**3/(Na*Nb))*histgofrAB[grbin]/(shellvol*meas)
    rmid = rinner + 0.5*grdelta
    print(rmid,gofrAA,gofrBB,gofrAB,file=fileoutgofr)
    
# plot T
tarray = np.arange(dt,(nstepMD+1)*dt,dt)
plt.figure()
fax = plt.axes()
fax.set_xlim(0,(nstepMD+1)*dt)
fax.set_xticks(np.arange(0,(nstepMD+1)*dt,nstepBoltz*dt))
fax.plot(tarray,temp_array,color='black')
fax.scatter(tarray,temp_array,s=40,color='blue')
fax.set_xlabel('$t$',fontsize=20)
fax.set_ylabel('$T_\mathrm{meas}$',fontsize=20)
fax.xaxis.set_major_locator(MultipleLocator(nstepBoltz*dt))
fax.xaxis.set_minor_locator(MultipleLocator(dt))
# plt.savefig('Toft_trial.eps') replace plt.show with this line
plt.show()

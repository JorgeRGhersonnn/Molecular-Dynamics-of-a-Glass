#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

# Parameters 
L = 9.4  # Box length
Ldiv2 = L / 2.0
dr = 0.1  # Bin width for RDF
Na = 80
Nb = 20
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

# Load data
data = np.genfromtxt("trajectory_data.txt", delimiter=",", skip_header=1, dtype=float, 
                     unpack=True)

unique_times = np.unique(data[:, 0])
n_times = len(unique_times)

if data.shape[0] % N != 0: # Check number of rows is a multiple of N
    raise ValueError("The number of rows in the data file is not a multiple of N.")

data_reshaped = data.reshape((n_times, N, -1))  # Reshape data to (n_times, N, -1)

# Extract x, y, z positions
positionsxyz = data_reshaped[:,:, 2:5]

# radial distribution function (RDF) analysis
def compute_rdf(positions, L, dr):
    """
    Compute the radial distribution function g(r) from positions.

    Parameters:
        positions (ndarray): Array of shape (N, 3) with x,y,z positions.
        L (float): Box length.
        dr (float): Bin width.

    Returns:
        r_vals (ndarray): Midpoints of the r bins.
        g_r (ndarray): Computed g(r) values.
    """
    N = positions.shape[0]
    dists = []
    # Loop over all unique pairs of particles
    for i in range(N-1):
        for j in range(i+1, N):
            dx = positions[j, 0] - positions[i, 0]
            dy = positions[j, 1] - positions[i, 1]
            dz = positions[j, 2] - positions[i, 2]
            # Apply minimum image convention:
            dx -= L * np.rint(dx / L)
            dy -= L * np.rint(dy / L)
            dz -= L * np.rint(dz / L)
            r = np.sqrt(dx*dx + dy*dy + dz*dz)
            dists.append(r)
    dists = np.array(dists)
    
    # Set up histogram bins
    bin_edges = np.arange(0.0, Ldiv2 + dr, dr)
    hist, bin_edges = np.histogram(dists, bins=bin_edges)
    r_vals = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    # Shell volumes
    shell_vol = (4.0/3.0) * np.pi * (np.power(bin_edges[1:], 3) - np.power(bin_edges[:-1], 3))
   
    # Normalize: ideal gas prediction (number of pairs expected in shell)
    n_pairs = N*(N-1)/2
    shell_vol = (4.0/3.0) * np.pi * (np.power(bin_edges[1:], 3) - np.power(bin_edges[:-1], 3))
    n_ideal = n_pairs * shell_vol / (L**3)
    # Avoid division by zero; note that very small r bins are often empty
    g_r = hist / n_ideal
    return r_vals, g_r

positions_final = positionsxyz[-1]  # Last snapshot
dr = 0.1
r_vals, g_r = compute_rdf(positions_final, L, dr)

# Figure for RDF

plt.figure(figsize=(7,5))
plt.plot(r_vals, g_r, 'bo-', label=f"t = {unique_times[-1]:.3f}")
plt.xlabel("r")
plt.ylabel("g(r)")
plt.title("Radial Distribution Function")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Initial positions
positions0 = positionsxyz[0]  # First snapshot

# Mean square displacement (MSD) analysis
msd=[]
for i in range(N):
    pos_t = positionsxyz[i]
    
    disp = pos_t - positions0
    disp_sq = np.sum(disp**2, axis=1)
    msd.append(np.mean(disp_sq))
    
msd = np.array(msd)

# Plot MSD
plt.figure(figsize=(7,5))
plt.plot(unique_times, msd, 'ro-', label="MSD")
plt.xlabel("Time")
plt.ylabel("MSD")
plt.title("Mean Square Displacement vs. Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

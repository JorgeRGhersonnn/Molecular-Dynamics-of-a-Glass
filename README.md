# MDGlass

## Equations of Motion

The motion of particles in a molecular system is governed by Newton’s second law, which states that the force acting on a particle is equal to its mass multiplied by its acceleration. In mathematical terms, for a system of $N$ particles, each having mass $m_i$, and position $\mathbf{r}_i\$, the equation takes the form:

$$\mathbf{F}_i = m_i \frac{d^2 \mathbf{r}_i}{dt^2}.$$

The force $\mathbf{F}_i$ is derived from the interatomic potential $U(\mathbf{r})$ by taking its negative gradient:

$$\mathbf{F}_i = -\mathbf{\nabla}_i \, U(\mathbf{r}_i.)$$

## The Lennard-Jones Potential
A fundamental model for describing interatomic interactions is the Lennard-Jones potential, which captures the balance between attractive and repulsive forces. This potential is mathematically represented as:

$$U(r) = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right],$$ 

where $\epsilon$ represents the depth of the potential well (i.e., the strength of the attractive interaction), $\sigma$ is approximately the distance at which the interparticle potential becomes zero, the $r^{-12}$ term models the steep repulsion due to overlapping electron orbitals (Pauli exclusion principle), and finally, the $r^{-6}$ term captures the long-range van der Waals attractions.

# Molecular Dynamics of a Glass

## Equations of Motion

The motion of particles in a molecular system is governed by Newton’s second law, which states that the force acting on a particle is equal to its mass multiplied by its acceleration. In mathematical terms, for a system of $N$ particles, each having mass $m_i$, and position $\mathbf{r}_i\$, the equation takes the form:

$$\mathbf{F}_i = m_i \frac{d^2 \mathbf{r}_i}{dt^2}.$$

The force $\mathbf{F}_i$ is derived from the interatomic potential $U(\mathbf{r})$ by taking its negative gradient:

$$\mathbf{F}_i = -\mathbf{\nabla}_i \, U(\mathbf{r}_i.)$$

## The Lennard-Jones Potential
A fundamental model for describing interatomic interactions is the Lennard-Jones potential, which captures the balance between attractive and repulsive forces. This potential is mathematically represented as:

$$U(r)=4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right],$$ 

where $\epsilon$ represents the depth of the potential well (i.e., the strength of the attractive interaction), $\sigma$ is approximately the distance at which the interparticle potential becomes zero, the $r^{-12}$ term models the steep repulsion due to overlapping electron orbitals (Pauli exclusion principle), and finally, the $r^{-6}$ term captures the long-range van der Waals attractions.

## Kob-Andersen Potential
The Kob-Andersen potential modifies the Lennard-Jones potential by allowing the parameters $\epsilon$ and $\sigma$ to depend on the particle types $\alpha,\beta \in \{A,B\}$ of  particles $i$ and $j$ for an $80:20 \, A:B$ mixture. This adjustment can be used to capture additional anharmonic effects in atomic vibrations. The potential is given by:

$$U_{\alpha\beta}(r_{ij}) = 4\epsilon_{\alpha\beta} \left[ \left(\frac{\sigma_{\alpha\beta}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{\alpha\beta}}{r_{ij}}\right)^{6}\right]$$

with the understanding that in some implementations, an additional Morse potential term may be included to further capture anharmonic effects. Here, we use the same units as in (Vollmayr-Lee et al., 2020), such that $\sigma_{AA} = 1$ (length unit), $\epsilon_{AA}=1$ (energy unit), $m_A = 1$ (mass unit), and $k_B=1$ for the temperature unit $\epsilon_{AA}/k_b$. The time unit is $\sqrt{m_A \sigma^2_{AA}/\epsilon_{AA}}$. This results in our KA parameters being $\sigma_{AA}=1,\, \epsilon_{AA}=1,\, \sigma_{AB}=0.8,\, \epsilon_{AB}=1.5,\, \sigma_{BB}=0.88,\, \epsilon_{BB}=0.5$, and $m_A = m_B = 1.0$.

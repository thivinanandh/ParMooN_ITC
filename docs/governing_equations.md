# Governing Equations

The complete model consists of three coupled systems: fluid flow equations, particle transport equation, and particle drift velocity equation. Each system is described in detail below.

# Governing Equations

## Non-dimensionalization

The governing equations are presented in non-dimensional form. The following non-dimensional quantities are introduced (denoted with tilde):

$$\tilde{x} = \frac{x}{L}, \quad \tilde{\mathbf{u}} = \frac{\mathbf{u}}{U}, \quad \tilde{\mathbf{u}}_b = \frac{\mathbf{u}_b}{U}, \quad \tilde{\mathbf{w}} = \frac{\mathbf{w}}{U}, \quad \tilde{t} = \frac{tU}{L}, \quad \tilde{p} = \frac{p}{\rho U^2}$$

where:
- $L$ is the characteristic length
- $U$ is the characteristic velocity
- Variables with tilde ($\tilde{}$) represent non-dimensional quantities

The relevant non-dimensional numbers are:

$$\text{Re} = \frac{\rho UL}{\mu_0} \quad \text{(Reynolds number)}$$

$$\text{We} = \frac{\rho U^2L}{\sigma} \quad \text{(Weber number)}$$

$$\text{Fr} = \frac{U^2}{Lg} \quad \text{(Froude number)}$$

$$\beta = \frac{1}{\varepsilon_{\mu}\rho U} \quad \text{(Slip parameter)}$$



## 1. Fluid Flow Equations

The incompressible fluid flow in the three-dimensional domain Ω is governed by the time-dependent Navier-Stokes equations:

$$\frac{\partial \mathbf{u}}{\partial t} - \frac{2}{\text{Re}}\nabla\cdot\mathbb{D}(\mathbf{u}) + (\mathbf{u}\cdot\nabla)\mathbf{u} + \nabla p = \mathbf{0}$$
$$\nabla\cdot\mathbf{u} = 0$$
$$\mathbf{u}(0,\mathbf{x}) = \mathbf{u}_0$$

where:
- $\mathbf{u}$ is the fluid velocity
- $p$ is the pressure
- $\text{Re}$ is the Reynolds number
- $\mathbb{D}(\mathbf{u})$ is the velocity deformation tensor
- $T$ is the final time
- $\mathbf{u}_0$ is the initial velocity

### Boundary Conditions
The system is completed with the following boundary conditions:

$$\mathbf{u}(t,\mathbf{x}) = \mathbf{u}_D \quad \text{on } (0,T] \times \Gamma_\text{in}$$
$$\mathbf{u}(t,\mathbf{x}) = \mathbf{0} \quad \text{on } (0,T] \times \Gamma_\text{wall}$$
$$\left(\frac{2}{\text{Re}}\mathbb{D}(\mathbf{u}) - p\mathbb{I}\right)\cdot\mathbf{n} = \mathbf{0} \quad \text{on } (0,T] \times \Gamma_\text{out}$$

where:
- $\mathbf{u}_D$ is the prescribed inlet velocity
- $\mathbb{I}$ is the identity tensor
- $\mathbf{n}$ is the outward normal to the boundary

The velocity deformation tensor and Reynolds number are defined as:
$$\mathbb{D}(\mathbf{u}) = \frac{\nabla\mathbf{u} + \nabla\mathbf{u}^T}{2}$$
$$\text{Re} = \frac{\rho UL}{\mu}$$

## 2. Population Balance Equation

The evolution and transport of particles is governed by a population balance equation (PBE) that accounts for both physical space coordinates and internal coordinates (representing particle properties such as size). The complete equation is:

$$\frac{\partial c}{\partial t} - \varepsilon\Delta_x c + \mathbf{u}_p \cdot \nabla_x c + \mathbf{g} \cdot \nabla_\ell c = f \quad \text{in } (0,T] \times \Omega$$

### Domain Definition

The computational domain $\Omega$ is defined as a Cartesian product of physical and internal coordinate spaces:

$$\Omega := \Omega_X \times \Omega_L \subset \mathbb{R}^d \times \mathbb{R}^s$$

where:
- $\Omega_X \subset \mathbb{R}^d$ is the physical space domain ($d = 2$ or $3$)
- $\Omega_L \subset \mathbb{R}^s$ is the internal coordinate domain ($s \geq 1$)
- $\partial\Omega$ represents the polyhedral boundary

### Equation Components

1. Temporal Term:
   $$\frac{\partial c}{\partial t}$$

2. Diffusion in Physical Space:
   $$\varepsilon\Delta_x c$$
   where $\Delta_x$ is the Laplace operator in physical space

3. Convection in Physical Space:
   $$\mathbf{u}_p \cdot \nabla_x c$$
   where $\nabla_x$ is the gradient operator in $\Omega_X$

4. Internal Coordinate Transport:
   $$\mathbf{g} \cdot \nabla_\ell c$$
   where $\nabla_\ell$ is the gradient operator in $\Omega_L$

### Boundary and Initial Conditions

$$c(t,x,\ell) = 0 \quad \text{on } (0,T] \times \partial\Omega$$
$$c(0,x,\ell) = c_0(x,\ell) \quad \text{in } \Omega$$

The internal coordinate $\ell$ represents particle properties (typically size or mass), making this a multi-dimensional problem that captures both spatial transport and the evolution of particle size distribution.

## 3. Particle Drift Velocity Equation

The particle drift velocity $\mathbf{u}_p$ is governed by the equation of motion:

$$m_p\frac{d\mathbf{u}_p}{dt} = \mathbf{F}_D + \mathbf{F}_G + \mathbf{F}_B$$

### Force Components

1. Drag Force:
   $$\mathbf{F}_D = \frac{3}{4}\frac{\rho_f}{\rho_p}\frac{m_p}{d_p}\frac{C_D}{C_C}|\mathbf{u}_f - \mathbf{u}_p|(\mathbf{u}_f - \mathbf{u}_p)$$

2. Gravitational Force:
   $$\mathbf{F}_G = m_p\mathbf{g}\frac{\rho_p - \rho_f}{\rho_p}$$

3. Brownian Force:
   $$\mathbf{F}_B \text{ (stochastic force term)}$$

where:
- $\rho_f$, $\rho_p$ are fluid and particle densities
- $d_p$ is particle diameter
- $C_D$ is drag coefficient
- $C_C$ is Cunningham correction factor
- $\mathbf{g}$ is gravitational acceleration
- $\mathbf{u}_f$ is fluid velocity

This system of equations is coupled through:
- The fluid velocity field affects particle transport through drift velocity
- The particle concentration can influence fluid properties
- The drift velocity depends on both fluid velocity and particle properties
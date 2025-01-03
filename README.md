# Eulerian Paricle Deposition for Human Air Pathways
---


A high-performance C++ implementation for simulating 3D particle deposition using an Eulerian framework. The solver couples incompressible Navier-Stokes equations with population balance modeling using finite element methods.

## Overview

This code provides a comprehensive solution for particle deposition problems by integrating:
- 3D incompressible fluid flow simulation
- Population Balance Equation (PBE) modeling
- Intel Pardiso solver integration for efficient computation
- Finite element methodology with operator splitting approach

## Key Features

- Fully coupled Eulerian framework for particle-fluid interaction
- Advanced finite element implementation for spatial discretization
- Intel Pardiso solver integration for enhanced performance
- Parallel computing capabilities
- Comprehensive post-processing tools

## Documentation

### [Governing Equations](docs/governing_equations.md)
- Fluid Flow Equations
  - Incompressible Navier-Stokes
  - Turbulence modeling
- Population Balance Equations
  - Transport equations
  - Source terms
  - Boundary conditions

### Test Cases

#### [Bent Pipe Case](docs/bent_pipe/README.md)
- Problem Description
  - Geometry and Setup
  - Boundary Conditions
  - Initial Conditions
- Results
  - Fluid Flow Results
    - Velocity profiles
    - Pressure distributions
  - PBE Results
    - Particle concentration
    - Deposition patterns
  - Validation with Literature

#### [Actual Geometry](docs/application/README.md)
- Problem Description
  - Geometry and Setup
  - Boundary Conditions
  - Initial Conditions
- Results
  - Fluid Flow Results
    - Velocity profiles
    - Pressure distributions
  - PBE Results
    - Particle concentration
    - Deposition patterns
  - Validation with Experiments

### [Code Implementation](docs/code_workflow.md)
- Main algorithm structure
- Class hierarchies
- Solver implementation


### [Installation and Setup](docs/installation.md)
- System requirements
- Dependencies
  - Intel Pardiso Solver
  - Other libraries
- Build instructions
- Running test cases

## References

1. Wilbrandt, U., Bartsch, C., Ahmed, N., Alia, N., Anker, F., Blank, L., Caiazzo, A., Ganesan, S., Giere, S., Matthies, G., Meesala, R., Shamim, A., Venkatesan, J., John, V. (2017). ParMooN—A modernized program package based on mapped finite elements. Computers and Mathematics with Applications.

2. Ganesan, S. (2012). An operator-splitting Galerkin/SUPG finite element method for population balance equations: Stability and convergence. ESAIM: Mathematical Modelling and Numerical Analysis (M2AN), 46, 1447-1465.

3. Ganesan, S., Tobiska, L. (2012). An operator-splitting finite element method for the efficient parallel solution of multidimensional population balance systems. Chemical Engineering Science, 69(1), 59-68.

## Languages and Tools
- C++
- Intel Pardiso
- Cmake Tools

## Authors
 - Thivin Anandh
 - Sashikumaar Ganesan

## Contact
 - SashiKumaar Ganesan (sashi@iisc.ac.in)
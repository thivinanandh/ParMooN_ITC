# Bent Pipe Case Study

## Geometry Description

The computational domain consists of a three-dimensional circular pipe with the following specifications:
- Diameter (D) = 10 mm
- Bend angle = 90 degrees
- Inlet section length = D
- Outlet section length = 2D
- Bend radius = 2.8D

![Bent Pipe Geometry](geometry.png)

The flow enters horizontally from the inlet at the top-left and exits vertically at the bottom-right. The geometry was generated using FreeCAD [1] and meshed using GMSH [2].

![Reference Results](reference.png)
*Reference results from AeroSolved implementation [3]*

## Mesh Statistics
```
Number of root Vertices: 15,480
Number of root cells: 71,423
Number of Faces: 149,566
Total Number of Cells: 71,423
```

## Physical Parameters

### Flow Properties
The continuous phase properties are:
- Dynamic viscosity (μc) = 10⁻⁵ kg/(ms)
- Density (ρc) = 1 kg/m³

### Non-dimensional Parameters
Using characteristic scales:
- Length scale (L) = 0.01 m
- Velocity scale (U) = 1 m/s

The Reynolds number is calculated as:

$$\text{Re} = \frac{\rho U L}{\mu} = \frac{1 \cdot 1 \cdot 0.01}{10^{-5}} = 1000$$

## Numerical Method

### FEM Discretization
- Velocity field: P2 elements (second order)
- Pressure field: P1 elements (first order)
- Number of DOF for velocity field: 109,102
- Matrix properties:
  - Dimensions: 109,102 × 109,102
  - Number of entries: 2,359,517
- Quadrature: P5 formula for tetrahedra

### Boundary Conditions
1. Inlet (left boundary):
   - Uniform velocity: U = (1, 0, 0) m/s

2. Outlet (bottom boundary):
   - Zero Neumann condition

3. Walls:
   - No-slip condition: u = 0

## References

[1] FreeCAD. (Version X.X) [Computer software]. https://www.freecadweb.org/

[2] Geuzaine, C., & Remacle, J.-F. (2009). Gmsh: A 3-D finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering.

[3] Lucci, F., Frederix, E. M. A., & Kuczaj, A. K. (2022). AeroSolved: Computational fluid dynamics modeling of multispecies aerosol flows with sectional and moment methods. Journal of Aerosol Science, 159, 105854.
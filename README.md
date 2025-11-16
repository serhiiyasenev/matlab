# VisualizeEllipsoidPotential

`VisualizeEllipsoidPotential` is a MATLAB function for visualizing the gravitational potential over the surface of a triaxial ellipsoid. The function simulates the distribution of gravitational potential on the surface of a celestial body with specified semi-axes and mass using spherical coordinates.

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=serhiiyasenev/matlab&file=https://github.com/serhiiyasenev/matlab/blob/master/VisualizeEllipsoidPotential.m)

## Overview

This script performs the following tasks:
- Generates a 3D ellipsoid defined by semi-axes `a`, `b`, and `c`.
- Computes the gravitational potential on the ellipsoid surface using a quadrupole approximation.
- Visualizes the 3D surface of the ellipsoid, colored by the computed gravitational potential.
- Plots equipotential contours directly on the ellipsoid surface for enhanced interpretation.

## Parameters

The body parameters are defined within the script:

| Variable | Description                        | Example Value   |
|----------|------------------------------------|-----------------|
| `m`      | Mass of the body (kg)              | `1481.90e20`    |
| `a`      | Semi-axis along x (m)              | `2634.40e3`     |
| `b`      | Semi-axis along y (m)              | `2634.10e3`     |
| `c`      | Semi-axis along z (m)              | `2633.80e3`     |
| `G`      | Gravitational constant (m³/kg/s²)  | `6.673e-11`     |

## Potential Formula

The gravitational potential is computed using the following approximation:

Φ(r, θ, φ) = G·m / r · [1 + ((2a²−b²−c²)x² + (2b²−a²−c²)y² + (2c²−a²−b²)z²) / (10·r⁴)]

This expression accounts for deviations from spherical symmetry due to the ellipsoidal shape and incorporates second-order (quadrupole-like) corrections.

## Visualization

- The ellipsoid surface is rendered using `surf` and shaded according to potential values.
- Equipotential lines are overlaid using `contourc` and plotted in 3D space.
- The view is 3D with equal axis scaling and interactive rotation enabled.
- A colorbar provides a scale for interpreting the potential magnitude.

## Requirements

- MATLAB (recommended version R2018 or newer)
- No external toolboxes required — uses core MATLAB graphics functions

## How to Run

To execute the visualization, simply run the following command in the MATLAB command window:

```matlab
VisualizeEllipsoidPotential
```

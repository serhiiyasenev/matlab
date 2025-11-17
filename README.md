# VisualizeEllipsoidPotential

`VisualizeEllipsoidPotential` is a MATLAB function for visualizing the gravitational and effective potential over the surface of a triaxial ellipsoid. The function simulates the distribution of gravitational potential on the surface of a celestial body (Ganymede) with specified semi-axes, mass, and rotation using spherical coordinates.

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=serhiiyasenev/matlab&file=https://github.com/serhiiyasenev/matlab/blob/master/VisualizeEllipsoidPotential.m)

## Overview

This script performs the following tasks:
- Generates a 3D ellipsoid defined by semi-axes `a`, `b`, and `c` (representing Ganymede).
- Computes the gravitational potential on the ellipsoid surface using a quadrupole approximation.
- Includes centrifugal effects from the body's rotation to calculate the effective potential.
- Visualizes the 3D surface with:
  - Gray ellipsoid base showing the physical surface
  - Color overlay (blue→yellow gradient) representing effective potential values
  - Black dashed isolines showing gravitational potential contours
  - Red arrows indicating the direction of gravitational force (centripetal)
- Shows how rotation affects the potential: the equator appears "higher" in potential due to centrifugal force.

## Parameters

The body parameters are defined within the script:

| Variable | Description                        | Example Value   |
|----------|------------------------------------|-----------------|
| `m`      | Mass of the body (kg)              | `1481.90e20`    |
| `a`      | Semi-axis along x (m)              | `2634.40e3`     |
| `b`      | Semi-axis along y (m)              | `2634.10e3`     |
| `c`      | Semi-axis along z (m)              | `2633.80e3`     |
| `G`      | Gravitational constant (m³/kg/s²)  | `6.673e-11`     |
| `omega`  | Angular velocity (rad/s)           | `1.016e-5`      |

## Potential Formula

The effective potential is computed as the sum of gravitational and centrifugal potentials:

### Gravitational Potential (Quadrupole Approximation):

Φ_grav(r, θ, φ) = G·m / r · [1 + ((2a²−b²−c²)x² + (2b²−a²−c²)y² + (2c²−a²−b²)z²) / (10·r⁴)]

This expression accounts for deviations from spherical symmetry due to the ellipsoidal shape and incorporates second-order (quadrupole-like) corrections.

### Centrifugal Potential:

Φ_centrifugal = -0.5 · ω² · (x² + y²)

Where ω is the angular velocity of rotation, and (x² + y²) represents the squared distance from the rotation axis.

### Effective Potential:

Φ_eff = Φ_grav + Φ_centrifugal

The effective potential shows that the equator has a slightly less negative (higher) potential than the poles due to the outward centrifugal force from rotation.

## Visualization

The visualization includes multiple layers to provide a comprehensive view of Ganymede's gravitational field:

- **Gray Ellipsoid Surface**: The base layer shows the physical ellipsoid surface in gray (semi-transparent).
- **Effective Potential Overlay**: A colored surface layer using a blue→yellow gradient represents the effective potential values. Blue indicates more negative (deeper) potential, while yellow shows less negative (higher) potential.
- **Black Dashed Isolines**: Equipotential lines for the gravitational potential are overlaid on the surface. These lines are nearly horizontal, converging toward the poles where the potential is deeper.
- **Red Arrows**: Vectors showing the direction of gravitational force, all pointing centripetally toward the body's center. Arrow bases are positioned slightly above the surface, with tips on the surface for visual clarity.
- **Physical Interpretation**: The visualization clearly shows that the equator has a slightly higher potential (less negative) than the poles due to centrifugal effects from rotation.

The view is 3D with equal axis scaling and interactive rotation enabled. A colorbar provides a scale for interpreting the effective potential magnitude.

## Requirements

- MATLAB (recommended version R2018 or newer)
- No external toolboxes required — uses core MATLAB graphics functions

## How to Run

To execute the visualization, run the following command in the MATLAB command window:

```matlab
VisualizeEllipsoidPotential
```

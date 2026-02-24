# Lattice Boltzmann Method (LBM) Fluid Simulation

A high-performance Python implementation of the Lattice Boltzmann Method (LBM) to simulate 2D fluid flow around obstacles. This solver demonstrates how microscopic particle interactions can recover macroscopic fluid dynamics.

## Overview
Unlike traditional CFD solvers that discretize the Navier-Stokes equations, this project uses the **Lattice Boltzmann Method** with the **D2Q9** lattice configuration (2 Dimensions, 9 Velocity vectors). 

The simulation implements the **BGK (Bhatnagar-Gross-Krook)** collision operator to model fluid behavior and visualize complex flow patterns, such as vortex shedding, around various geometries.

## Key Features
* **D2Q9 Lattice Model:** Standard 2D grid with 9 discrete velocity directions.
* **Geometry Interaction:** Boundary conditions for flow around cylinders, squares, and arbitrary obstacles.
* **Vortex Shedding:** Capable of simulating the von Kármán vortex street at different Reynolds numbers.
* **Efficient NumPy Implementation:** Vectorized operations for faster computation of the streaming and collision steps.



## The Physics Behind It
The solver follows the standard LBM algorithm:
1.  **Streaming:** Moving probability density functions to neighboring nodes.
2.  **Collision:** Relaxing the distributions toward equilibrium using the BGK operator.
3.  **Boundary Handling:** Implementing "bounce-back" conditions for solid obstacles.

## Getting Started
### Prerequisites
* Python 3.x
* NumPy
* Matplotlib (for visualization)

### Execution
Run the main script to start the simulation and view the real-time velocity field:
```bash
python lbm_solver.py
```
## Credits & References
This implementation is based on the work of Philip Mocz (2020), Princeton University.
* **Author**: Philip Mocs ([@PMocz](https://github.com/pmocz/))
* **Original Tutorial**: [Create Your Own Lattice Boltzmann Simulation with Python](https://medium.com/swlh/create-your-own-lattice-boltzmann-simulation-with-python-8759e8b53b1c)

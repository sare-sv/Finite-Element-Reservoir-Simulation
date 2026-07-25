# Finite Element Reservoir Simulation

## Overview

This project presents a one-dimensional finite element simulation of steady-state single-phase flow in a heterogeneous oil reservoir.

The objective is to investigate pressure distribution in a reservoir with heterogeneous permeability and determine the minimum water injection rate required to maintain reservoir pressure above the bubble point pressure.

The model is developed using Python and implements both linear and quadratic finite element formulations.

---

## Problem Description

The reservoir model consists of a 1D heterogeneous domain with two permeability zones.

Reservoir characteristics:

- Length: 9000 ft
- Cross-sectional area: 100,000 ft²
- Porosity: 0.18
- Oil viscosity: 5 cp
- Boundary pressure: 4000 psi

Permeability distribution:

- Zone 1: 60 md
- Zone 2: 20 md

Well configuration:

- Producer 1: 1500 ft
- Water injector: 4500 ft
- Producer 2: 7500 ft

---

## Numerical Method

The governing flow equation is solved using the finite element method.

Implemented approaches:

- Linear finite elements
- Quadratic finite elements

The global stiffness matrix is assembled from local element matrices, and Dirichlet boundary conditions are applied at both reservoir boundaries.

---

## Injection Optimization

An iterative search algorithm is used to determine the minimum water injection rate required to keep reservoir pressure above the bubble point pressure.

The injection rate is increased until:
Minimum reservoir pressure ≥ Bubble point pressure

---

## Results

For the selected mesh:

- Minimum required injection rate:
  - 1397 STB/day

- Minimum reservoir pressure after injection:
  - 2200.2 psi

- Material balance:
  - Total inflow = 2000 STB/day
  - Total outflow = 2000 STB/day

The linear and quadratic FEM solutions show very close agreement.

---

## Pressure Distribution

The following figure shows the pressure distribution obtained using linear and quadratic finite element formulations.

![Pressure Distribution](results/reservoir_pressure_plot.png)

---

## Requirements

Python libraries:

- NumPy
- Matplotlib

Install requirements using:
pip install numpy matplotlib

---

## Author

Sare Vahedi

M.Sc. Petroleum Engineering  
Reservoir Engineering

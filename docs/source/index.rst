.. pkmodel documentation master file, created by
   sphinx-quickstart on Fri Oct 31 12:11:20 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pkmodel documentation
=====================

pkmodel is a Python library for creating, solving, and visualising pharmacokinetic (PK) models.
It allows the user to represent an organism as a set of interacting compartments and simulate the movement of substances between them over time. 

Features
==========
Flexible compartment modelling 

- Define custom compartments to represent physiological spaces.

Mass transfer and elimination

- Add **fluxes** to describe one-way or bidirectional flows between compartments. 
- Add **clearances** to describe the elimination process.
- Define **dosage regimens** to represent substance administration (continuous, bolus, or custom time-dependent inputs)

Rate laws

- Both **first-order** and **zero-order** kinetics supported, in addition to custom time-dependent dosage regimens.

Automatic ODE construction

- Builds a system of linear (possibly time-dependent) ordinary differential equations that describe the model dynamics.

Simulation and visualisation

- Provides built-in numerical solvers and visualisation tools for plotting the results of the simulations.
- Also includes tools that can generate a graphical overview of model architecture and a tool to generate a Markdown representation of the ODE system. 


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   pkmodel
   pkmodel.tests


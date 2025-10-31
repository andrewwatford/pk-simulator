![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Build](https://img.shields.io/github/actions/workflow/status/yourusername/pk-simulator/tests.yml)
![Python](https://img.shields.io/badge/python-3.9%2B-blue)
# pkmodel: Pharmacokinetic Modelling Library

pkmodel is a Python library for creating, solving, and visualising pharmacokinetic (PK) models.
It allows the user to represent an organism as a set of interacting compartments and simulate the movement of substances between them over time. 

## Table of Contents
- [Library Features](#library-features)
- [Installation](#installation)
- [Usage](#usage)
- [Configuration Example](#alternatively-the-model-can-be-instantiated-using-a-config-dictionary)
- [License](#license)
- [Contributions](#contributions)

## Library Features
### Flexible compartment modelling 
- Define custom compartments to represent physiological spaces.
### Mass transfer and elimination
- Add **fluxes** to describe one-way or bidirectional flows between compartments. 
- Add **clearances** to describe the elimination process.
- Define **dosage regimens** to represent substance administration (continuous, bolus, or custom time-dependent inputs)
### Rate laws
- Both **first-order** and **zero-order** kinetics supported, in addition to custom time-dependent dosage regimens.
### Automatic ODE construction
- Builds a system of linear (possibly time-dependent) ordinary differential equations that describe the model dynamics.
### Simulation and visualisation
- Provides built-in numerical solvers and visualisation tools for plotting the results of the simulations.
- Also includes tools that can generate a graphical overview of model architecture and a tool to generate a Markdown representation of the ODE system. 


## Installation
### via pip
     pip install pkmodel

### via source
     git clone https://github.com/andrewwatford/pk-simulator.git
     cd pk-simulator
     pip install –e

## Usage

### Instantiate the model
```python
import pkmodel as pk

model = pk.CompartmentModel()
```

### Add compartments to the model
```python
central = Compartment(volume=22)
peripheral = Compartment(volume=7)

model.add_compartment(central)
model.add_compartment(peripheral)
```

### Add fluxes
```python
c_p_flux = Flux(
        source=central,
        dest=peripheral,
        nature="bidirectional",
        rate_law="first",
        rate_constant=5,
    )

model.add_flux(c_p_flux)
```

### Add dosages and clearances
```python
central_clr = Clearance(
        source=central,
        rate_constant=5,
        rate_law="first"
    )

central_dsg = Dosage(
    dest=central,
    regime="constant",
    rate_constant=1
)

model.add_dosage(central_dsg)
model.add_clearance(central_clr)
```
### Run the model and plot the results
```python
# Define initial conditions and volumes
y0 = [0, 0]  # Initial mass in each compartment
t_span = [0, 30]  # Time span for the simulation

# Run the simulation
result = model.run(t_span, y0)

# Plots
import matplotlib.pyplot as plt

fig, axs = model.plot_all(result)
plt.savefig('./example.png')
```
## Visulaization
### Via Graphviz 
Note: Graphviz needs to be installed from your system's package manager first
```python
pip install graphviz
Fig = model.plot_using_graphviz(filename="compartmentmodel_graphviz")
plt.savefig("compartmentmodel_graphviz.png")
```

## Alternatively, the model can be instantiated using a config dictionary
```python
import pkmodel as pk
import matplotlib.pyplot as plt


# Define model config
config = {

    "compartments": {
        "central": 22.0,
        "peripheral": 7.0,
    },

    "fluxes": {
        "c_p": {
            "source":"central",
            "dest": "peripheral",
            "rate_constant": 5.0,
            "nature":"bidirectional",
            "rate_law":"first"
        }
    },

    "clearances": {
        "central_clearance":{
            "source":"central",
            "rate_constant": 5.0,
            "rate_law":"first"
        }
    },

    "dosages": {
        "central_dosage":{
            "dest":"central",
            "regime":"constant",
            "rate_constant": 1.0,
        }
    }
}

#Instantiate the model object
model = pk.CompartmentModel.from_config(config)

#Build and solve the model
model.build_linear_rhs()

#Initial conditions and volumes
y0 = [0, 0]  # Initial mass in each compartment
t_span = [0, 30]  # Time span for the simulation

#Run the simulation
result = model.run(t_span, y0)

#Plots
fig, axs = model.plot_all(result)
plt.savefig('./example.png')
```

## License
This project is licensed under the MIT license
## Contributions

To contribute or fix an issue, please open an issue or submit a pull request on GitHub.

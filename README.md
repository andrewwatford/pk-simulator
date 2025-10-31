![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Build](https://img.shields.io/github/actions/workflow/status/andrewwatford/pk-simulator/.github/workflows/CI_workflow.yml)
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
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

## Alternatively, the model can be instantiated using a config file
```python
model4 = CompartmentModel.from_json("pkmodel/config.json")
```
Note that at present, time-dependent dosages are only able to be added by creating a `Dosage` object directly, and adding it to the model. For example:
```python
central_dsg = Dosage(
    id='central_dsg',
    dest=central,
    regime='custom',
    dosage_func=my_dosage_func
)
model.add_dosage(central_dsg)
```
Also note that editing the attributes of a `CompartmentModel` instance or its component instances should be done at the user's own risk and may lead to unexpected behaviour (e.g., fluxes that are not linked to existing compartments).

## Visualisation
### Printing out the model architecture
```python
print(model)
```
```text
This is a model containing the following compartments:
        Compartment 'central', with a volume of 22.0 L
        Compartment 'peripheral', with a volume of 7.0 L

These are connected by the following fluxes (if any):
        Flux 'c_p' (bidirectional, first-order, with a rate constant of 5.0), connecting compartments 'central' and 'peripheral'.

With the following dosages(if any):
        Dosage 'central_dosage' (constant), representing administration to the copmartment 'central'.

And the following clearances (if any):
        Clearance 'central_clearance' (first-order, with a rate constant of 5.0), representing elimination from the compartment 'central'.
```
### Printing the system of ODEs associated with the model
Calling the following method writes the system of linear ordinary differential describing the model dynamics to a .md file:
```python
model.generate_markdown("filename")
```
An example output looks like this:

$$
\frac{d q_{0}}{d t} = D_{0}(t)- C_{0}\frac{q_{0}}{V_{0}}- k_{0,1}\left(\frac{q_{0}}{V_{0}} - \frac{q_{1}}{V_{1}}\right)
$$
$$
\frac{d q_{1}}{d t} = - k_{0,1}\left(\frac{q_{1}}{V_{1}} - \frac{q_{0}}{V_{0}}\right)
$$

The .md file also contains a table relating numerical indices to compartment names, and a table defining each variable.

### Printing out model graphically
We also have rudimentary support for printing out a graphical representation of the model in two different ways (pyplot and graphviz).
#### Via Pyplot
A NetworkX object can be created using the `construct_graph` method. The `draw_basic_graph_pyplot` method can be used to produce the graphical representation. For example:
```
fig, axs = model.draw_basic_graph_pyplot()
```
#### Via Graphviz 
Note: Graphviz needs to be installed from your system's package manager (apt-get, etc.) first
```python
pip install graphviz
Fig = model.plot_using_graphviz(filename="compartmentmodel_graphviz")
plt.savefig("compartmentmodel_graphviz.png")
```

## License
This project is licensed under the MIT license
## Contributions

To contribute or fix an issue, please open an issue or submit a pull request on GitHub.

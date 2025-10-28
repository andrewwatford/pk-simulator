# pk-model Python Library

pk-model is a Python library for creating, solving, and visualising the solution of pharmokinetic models. 
## Features
- Specify compartments and volumes
- Add fluxes between compartments
- Add clearance between compartments
- Add dosages to compartments
- Specify initial conditions and solve the model over a specified time span
- Plot the results

## Installation
### The latest release can be installed via pip
     pip install pk-model

## Example analysis script:
---
     from pkmodel import CompartmentModel
     model = CompartmentModel(['Central', 'Peripheral'], [3.0, 5.0])
     model.add_flux('Central', 'Peripheral', rate_constant=0.5, rate_law='first')
     model.add_clearance('Central', rate_constant=0.3, rate_law='first')
     model.add_dosage('Central', lambda t: 10 if t < 1 else 0)
     t_span = (0, 10)
     y0 = [0, 0]
     sol = model.build(t_span, y0, t_eval=np.linspace(0, 10, 100)) 
     print(sol.t)
     print(sol.y)




# pkmodel Python Library

pk-model is a Python library for creating, solving, and visualising the solution of pharmokinetic models. 
## Library Features
- Specify compartments and volumes
- Add fluxes between compartments(bidirectional/unidirectional)
- Add clearance between compartments
- Specify rate law for fluxes and clearances
- Add dosages to compartments
- Build RHS of the ODE system
- Run simulations 
- Plot the results

## Installation
### via pip
     pip install pkmodel
### via source
     git clone https://github.com/andrewwatford/pk-simulator.git
     cd pk-simulator
     pip install –e

## Basic Example Script:
---
     import matplotlib.pyplot as plt
     from pkmodel.CompartmentModel import CompartmentModel
     
     #Define model config
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
     model = CompartmentModel.from_config(config)
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

## License
---
This project is licensed under the MIT license
## Contributions
---
To contibute, please open an issue or submit a pull request on GitHub.

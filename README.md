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
### The pkmodel library can be installed via pip
     pip install pkmodel

## Basic Example Script:
---
     * Needs to be updated with a better example
     #create compartments
     central = Compartment(id="central", volume=22)
     peripheral = Compartment(id="peripheral", volume=7)
     
     #create flux
     c_p_flux = Flux(id="c_p_flux", source=central, dest=peripheral, rate_constant=5, nature="bidirectional", rate_law="first")
     
     #create clearance
     central_clr = Clearance(id="central_clearance", source=central, rate_constant=5, rate_law="first")
     
     #create dosage
     central_dsg = Dosage(id="central_dosage", dest=central, regime="constant", rate_constant=1)
     
     #create model and add components
     model = CompartmentModel()
     model.add_compartment(central)
     model.add_compartment(peripheral) 
     model.add_flux(c_p_flux)
     model.add_clearance(central_clr)
     model.add_dosage(central_dsg)
     
     #build rhs
     model.build_linear_rhs()
     
     #run simulation  
     ds = model.run(t_span=(0, 10), y0=[22, 7], t_eval=np.linspace(0, 10, 100))
     
     #plot results
     fig, axs = model.plot_all(ds) 




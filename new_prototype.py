from pkmodel.builtin_fluxes import *

volumes_dict = {'central': 22, 'peripheral': 7}
model = CompartmentModel(compartment_names = list(volumes_dict.keys()), compartment_volumes = list(volumes_dict.values()))
model.add_flux(from_compartment = 'central', 
               to_compartment = 'peripheral', 
               rate_function = first_order_flux)
model.add_dosage(compartment = 'central', dosage_function = constant_dose)
model.add_clearance(compartment = 'central', clearance_function = first_order_clearance)
out = model.run(initial_conditions = [0, 0], final_time = 10)


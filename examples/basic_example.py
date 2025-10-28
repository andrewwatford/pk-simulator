import matplotlib.pyplot as plt
from pkmodel.builtin_fluxes import constant_dose
from pkmodel.CompartmentModel import CompartmentModel

# Define model
volumes_dict = {'central': 22, 'peripheral': 7}
model = CompartmentModel(compartment_names = list(volumes_dict.keys()), compartment_volumes = list(volumes_dict.values()))
# Adding fluxes, dosage, and clearances
model.add_flux(from_compartment = 'central', 
               to_compartment = 'peripheral', 
               rate_law = 'first', # or 'zero'
                rate_constant = 5)
model.add_dosage(compartment_name = 'central', dosage_func = constant_dose(1))
model.add_clearance(from_compartment = 'central', 
                    rate_law = 'first', # or 'zero'
                    rate_constant = 5)
# running the model
out = model.run(y0 = [0, 0], t_span = [0, 30])
# Plots
fig, axs = model.plot_all(out)
plt.savefig('./example.png')
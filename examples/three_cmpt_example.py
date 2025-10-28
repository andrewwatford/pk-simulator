import matplotlib.pyplot as plt
import numpy as np
from pkmodel.CompartmentModel import CompartmentModel

# Define model
volumes_dict = {'dose': 10, 'peripheral': 7, 'clearance': 20}
model = CompartmentModel(compartment_names = list(volumes_dict.keys()), compartment_volumes = list(volumes_dict.values()))
# Adding fluxes, dosage, and clearances
model.add_flux(from_compartment = 'dose',
                to_compartment = 'peripheral', 
                rate_law = 'first', # or 'zero'
                rate_constant = 3)
model.add_flux(from_compartment = 'peripheral',
                to_compartment = 'clearance',
                rate_law = 'first', # or 'zero'
                rate_constant = 2)
model.add_flux(from_compartment = 'dose',
                to_compartment = 'clearance',
                rate_law = 'first', # or 'zero'
                rate_constant = 2)
model.add_clearance(from_compartment = 'clearance', 
                    rate_law = 'first', # or 'zero'
                    rate_constant = 5)
# running the model
out = model.run(y0 = [1, 0, 0], t_span = [0, 30], t_eval = np.linspace(0, 30, 1000))
# Plots
fig, axs = model.plot_all(out)
plt.savefig('./three_cmpt_example.png')
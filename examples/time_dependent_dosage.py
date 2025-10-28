import matplotlib.pyplot as plt
import numpy as np
from pkmodel.CompartmentModel import CompartmentModel

def my_dosage(t):
    # custom dosage that spikes every three time units
    return 10 * np.cos(np.pi * t / 3)**10

# Define model
volumes_dict = {'central': 22, 'peripheral': 7}
model = CompartmentModel(compartment_names = list(volumes_dict.keys()), compartment_volumes = list(volumes_dict.values()))
# Adding fluxes, dosage, and clearances
model.add_flux(from_compartment = 'central', 
               to_compartment = 'peripheral', 
               rate_law = 'first', # or 'zero'
                rate_constant = 5)
model.add_dosage(compartment_name = 'central', dosage_func = my_dosage)
model.add_clearance(from_compartment = 'central', 
                    rate_law = 'first', # or 'zero'
                    rate_constant = 5)
# running the model
out = model.run(y0 = [0, 0], t_span = [0, 30], t_eval = np.linspace(0, 30, 1000))
# Plots
fig, axs = model.plot_all(out)
plt.savefig('./time_dependent_dosage.png')
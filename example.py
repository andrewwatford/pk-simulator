import matplotlib.pyplot as plt
from pkmodel.builtin_fluxes import constant_dose
from pkmodel.CompartmentModel import CompartmentModel

volumes_dict = {'central': 22, 'peripheral': 7}
model = CompartmentModel(compartment_names = list(volumes_dict.keys()), compartment_volumes = list(volumes_dict.values()))
model.add_flux(from_compartment = 'central', 
               to_compartment = 'peripheral', 
               rate_law = 'first', # or 'zero'
                rate_constant = 5)
model.add_dosage(compartment_name = 'central', dosage_func = constant_dose(1))
model.add_clearance(from_compartment = 'central', 
                    rate_law = 'first', # or 'zero'
                    rate_constant = 5)
# running the model
out = model.run(y0 = [50, 0], t_span = [0, 30])
# Plots
fig, axs = plt.subplots(2, 1, sharex=True)
out.central.plot(ax=axs[0])
axs[0].set_xlabel(None)
axs[0].set_ylabel('$q_C$ (central)')
out.peripheral.plot(ax=axs[1])
axs[1].set_xlabel('$t$')
axs[1].set_ylabel('$q_P$ (peripheral)')
axs[0].set_title('Compartment masses over time')
plt.savefig('./example.png')
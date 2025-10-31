import matplotlib.pyplot as plt
import numpy as np
from pkmodel.CompartmentModel import CompartmentModel, Compartment, Flux, Clearance, Dosage


def my_dosage(t):
    # custom dosage that spikes every three time units
    return 10 * np.cos(np.pi * t / 3)**10


# Define model compartments
central = Compartment(
    id='central',
    volume=22)
peripheral = Compartment(
    id='peripheral',
    volume=7)

# Define model flux, clearance, and dosage
c_p_flux = Flux(
    id='c_p_flux',
    source=central,
    dest=peripheral,
    rate_law='first',
    rate_constant=5,
    nature='bidirectional'
)

central_clr = Clearance(
    id='central_clr',
    source=central,
    rate_law='first',
    rate_constant=5
)

central_dsg = Dosage(
    id='central_dsg',
    dest=central,
    regime='custom',
    dosage_func=my_dosage
)

# Create compartment model
model = CompartmentModel()
model.add_compartment(central)
model.add_compartment(peripheral)
model.add_flux(c_p_flux)
model.add_clearance(central_clr)
model.add_dosage(central_dsg)

# Initial conditions, time span and evaluation points, and build and run
y0 = [0, 0]
t_span = [0, 30]
model.build_linear_rhs()
result = model.run(t_span, y0)

# Plots
fig, axs = model.plot_all(result)
plt.savefig('./time_dependent_dosage.png')
fig, axs = model.draw_basic_graph()
plt.savefig('./time_dependent_dosage_graph.png')

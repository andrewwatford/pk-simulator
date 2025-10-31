import matplotlib.pyplot as plt
from pkmodel.CompartmentModel import CompartmentModel

# Define model config
config = {

    "compartments": {
        "central": 22.0,
        "peripheral": 7.0,
    },

    "fluxes": {
        "c_p": {
            "source": "central",
            "dest": "peripheral",
            "rate_constant": 5.0,
            "nature": "bidirectional",
            "rate_law": "first"
        }
    },

    "clearances": {
        "central_clearance": {
            "source": "central",
            "rate_constant": 5.0,
            "rate_law": "first"
        }
    },

    "dosages": {
        "central_dosage": {
            "dest": "central",
            "regime": "constant",
            "rate_constant": 1.0,
        }
    }
}

# Instantiate the model object
model = CompartmentModel.from_config(config)
# Build and solve the model
model.build_linear_rhs()
# Initial conditions and volumes
y0 = [0, 0]  # Initial mass in each compartment
t_span = [0, 30]  # Time span for the simulation
# Run the simulation
result = model.run(t_span, y0)
# Plots
fig, axs = model.plot_all(result)
plt.savefig('./example.png')
fig, axs = model.draw_basic_graph()
plt.savefig('./example_graph.png')

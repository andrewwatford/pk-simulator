import matplotlib.pyplot as plt
from pkmodel.CompartmentModel import CompartmentModel

# Define model config
config = {

    "compartments": {
        "dose": 10.0,
        "peripheral": 7.0,
        "clearance": 20.0,
    },

    "fluxes": {
        "d_p": {
            "source": "dose",
            "dest": "peripheral",
            "rate_constant": 3.0,
            "nature": "bidirectional",
            "rate_law": "first"
        },
        "p_clr": {
            "source": "peripheral",
            "dest": "clearance",
            "rate_constant": 2.0,
            "nature": "bidirectional",
            "rate_law": "first"
        },
        "d_clr": {
            "source": "dose",
            "dest": "clearance",
            "rate_constant": 2.0,
            "nature": "bidirectional",
            "rate_law": "first"
        }
    },

    "clearances": {
        "clearance_clearance": {
            "source": "clearance",
            "rate_constant": 5.0,
            "rate_law": "first"
        }
    }
}
# Instantiate the model object
model = CompartmentModel.from_config(config)
# Build and solve the model
model.build_linear_rhs()
# Initial conditions and volumes
y0 = [1, 0, 0]
t_span = [0, 30]
# running the model
result = model.run(
    y0=y0,
    t_span=t_span)
# Plots
fig, axs = model.plot_all(result)
plt.savefig('./three_cmpt_example.png')
fig, axs = model.draw_basic_graph_pyplot()
plt.savefig('./three_cmpt_example_graph.png')

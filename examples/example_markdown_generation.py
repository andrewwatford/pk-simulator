import pkmodel as pk

config = {
        "compartments": {
            "central":    22.0,
            "peripheral": 7.0,
        },

        "fluxes": {
            "c_p": {
                "source":"peripheral",
                "dest": "central",
                "rate_constant": 25.0,
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
                "regime":"custom",
                "rate_constant": 0,
            }
        }
    }

model = pk.CompartmentModel.from_config(config)

model.generate_markdown("example_markdown")


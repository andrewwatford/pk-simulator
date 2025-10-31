import pkmodel as pk

# Example 1: The two-compartment intravenous bolus model given in the instructions
config_simple_intravenous_bolus = {
        "compartments": {
            "central":    22.0,
            "peripheral": 7.0,
        },

        "fluxes": {
            "c_p": {
                "source":"central",
                "dest": "peripheral",
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
                "rate_constant":0
            }
        }
    }

# Example 2: The three-compartment subcutaneous dosing model given in the instructions
config_simple_subcutaneous = {
        "compartments": {
            "absorbing": 5.0,
            "central":    22.0,
            "peripheral": 7.0,
        },

        "fluxes": {
            "a_c": {
                "source":"absorbing",
                "dest": "central",
                "rate_constant": 10,
                "nature":"unidirectional",
                "rate_law":"first"                
            },
            "c_p": {
                "source":"central",
                "dest": "peripheral",
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
            "absorbing_dosage":{
                "dest":"absorbing",
                "regime":"custom",
                "rate_constant":0
            }
        }
    }

# Example 3: An example of a physiologically-based model
# Central compartment is Arteries
# Gut, Liver, Pancreas, Heart/Lung, Brain, Muscle/Skin, Kidney have bidirectional diffusive exchange with arteries
# Liver & Gut have additional bidirectional diffusive exchange 
# Liver & Pancreas have additional bidirectional diffusive exchange
# Ingestion is via the Gut
config_physiological = {
        "compartments": {
            "ingestion":1,
            "gut":1,
            "arteries":1,
            "liver":1,
            "pancreas":1,
            "heart_and_lung":1,
            "brain":1,
            "muscle_and_skin":1,
            "kidney":1
        },

        "fluxes": {
            "ingestion_into_gut": {
                "source":"ingestion",
                "dest": "gut",
                "rate_constant": 1,
                "nature":"unidirectional",
                "rate_law":"first"
            },
            "liver_gut": {
                "source":"liver",
                "dest": "gut",
                "rate_constant": 1,
                "nature":"bidirectional",
                "rate_law":"first"
            },
            "liver_pancreas": {
                "source":"liver",
                "dest": "pancreas",
                "rate_constant": 1,
                "nature":"bidirectional",
                "rate_law":"first"
            },
        },

        "clearances": {
        },

        "dosages": {
            "ingestion_dosage":{
                "dest":"ingestion",
                "regime":"custom",
                "rate_constant":0
            }
        }
    }
# Add fluxes involving the central heart compartment
flux_dict = {}
for cname in config_physiological['compartments'].keys():
    if cname not in ['arteries','ingestion']:
        flux_dict[f'arteries_{cname}'] = {
                                            "source":"arteries",
                                            "dest": cname,
                                            "rate_constant": 1,
                                            "nature":"bidirectional",
                                            "rate_law":"first"
                                        }
config_physiological['fluxes'] = {**config_physiological['fluxes'], **flux_dict}


for name, config in zip(['intravenous_bolus', 'subcutaneous', 'physiologically_based'],
                        [config_simple_intravenous_bolus, config_simple_subcutaneous, config_physiological]):
    model = pk.CompartmentModel.from_config(config)
    model.generate_markdown(f"example_markdown_{name}")


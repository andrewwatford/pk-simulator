import pytest
import logging

import pkmodel as pk
from pkmodel.builtin_fluxes import constant_dose
import numpy as np

class TestCompartment:
    """
    Tests the CompartmentModel class.
    """
    def test_create(self):
        """
        Tests CompartmentModel creation.
        """
        # Instatiate a dict with compartment properties
        compartments_dict = {
            'central': 22, 
            'peripheral': 7}
        
        model = pk.CompartmentModel(
            compartment_names   = list(compartments_dict.keys()),
            compartment_volumes = list(compartments_dict.values()))
        
        # Check if attributes are stored correctly in a model object
        assert model.compartment_names == ["central", "peripheral"]
        assert model.compartment_volumes == [22, 7]

    @pytest.mark.parametrize(
    "name, flux_kwargs, clearance_kwargs, dose, y_init, expected" ,
        [
            (
                    "first_order_flux_clearance_dose",
                {
                    "from_compartment": "central",
                    "to_compartment": "peripheral",
                    "rate_constant": 5, 
                    "rate_law": "first"
                },
                {
                    "from_compartment": "central",
                    "rate_constant": 5,
                    "rate_law": "first"
                },
                1,  # Dose into the central compartment                    
                np.array([22.0, 7.0]),
                np.array([-4.0, 0.0])
            ),

            (
                    "first_order_flux_dose",
                {
                    "from_compartment": "central",
                    "to_compartment": "peripheral",
                    "rate_constant": 5, 
                    "rate_law": "first"
                },
                None,
                1,  # Dose into the central compartment                    
                np.array([22.0, 7.0]),
                np.array([1.0, 0.0])
            ),
        ]
    )

    def test_build_paramatrized(self, name, flux_kwargs, clearance_kwargs, dose, y_init, expected):
        logging.info(f"Testing the model {name}")
        model = pk.CompartmentModel(
            compartment_names   = ["central", "peripheral"],
            compartment_volumes = [22,7]
        )

        # Add optional flux
        if flux_kwargs is not None:
            model.add_flux(**flux_kwargs)

        # Add optional clearance
        if clearance_kwargs is not None:
            model.add_clearance(**clearance_kwargs)

        # Add dosage to central
        model.add_dosage(compartment_name="central", dosage_func=constant_dose(dose))

        # Build the model
        model.build()

        assert model.model_built == True
        assert model.model_changed_since_last_build == False

        # Testing if the RHS works with some sample numbers
        result = model.rhs(5, np.array(y_init))
        assert result == pytest.approx(np.asarray(expected), abs=1e-8)

                
                
            
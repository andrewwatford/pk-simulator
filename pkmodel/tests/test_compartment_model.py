import numpy.testing as npt
import pytest
import logging
import pkmodel as pk
from pkmodel.builtin_fluxes import constant_dose
import numpy as np

class TestCompartmentModel:
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

    def test_create_with_invalid_inputs(self):
        """
        Tests CompartmentModel creation with mismatched lengths of names and volumes.
        """

        # Test with mismatched lengths of names and volumes
        with pytest.raises(ValueError):
            pk.CompartmentModel(
                compartment_names   = ['central','peripheral'],
                compartment_volumes = [1,2,3])
    
    @pytest.mark.parametrize(
        "compartment_names, compartment_volumes, flux_dict_list, clearance_dict_list, expected_matrix, expected_cst_vector",
        [
            # Test case 1: two compartments, one first-order flux, no clearances
            (   
                ['central', 'peripheral'],
                [22, 7],
                [{
                    'from_compartment': 'central',
                    'to_compartment': 'peripheral',
                    'rate_constant': 1,
                    'rate_law': 'first'}],
                [{
                    'from_compartment': 'central',
                    'rate_constant': 0,
                    'rate_law': 'zero'
                }],
                [
                [-1/22, 1/7],
                [1/22, -1/7]
                ],
                [0, 0]
            ),
            # Test case 2: two compartments, one first-order flux, one first-order clearance,
            (   
                ['central', 'peripheral'],
                [22, 7],
                [{
                    'from_compartment': 'central',
                    'to_compartment': 'peripheral',
                    'rate_constant': 1,
                    'rate_law': 'first'},
                    ],
                [{
                    'from_compartment': 'central',
                    'rate_constant': 2,
                    'rate_law': 'first'
                }],
                [
                [-3/22, 1/7],
                [1/22, -1/7]
                ],
                [0, 0]
            ),
            # Test case 3: subcutaneous dosing. three compartments, 
            # Try to make edge cases for 2d and 3d
            # Make tests for 3d
        ]
    )

    def test_matrix_and_constant_vector(self, compartment_names, compartment_volumes, flux_dict_list, clearance_dict_list, expected_matrix, expected_cst_vector):
        """
        Tests that the RHS matrix and constant vector are correctly constructed
        after various combinations of zero- or first-order fluxes and clearances are added
        """      

        # Initialise compartment model
        model = pk.CompartmentModel(
            compartment_names   = compartment_names,
            compartment_volumes = compartment_volumes)

        # Add fluxes
        for flux_dict in flux_dict_list:
            model.add_flux(**flux_dict)
        
        # Add clearances
        for clearance_dict in clearance_dict_list:
            model.add_clearance(**clearance_dict)

        npt.assert_array_almost_equal(model.rhs_matrix, expected_matrix)
        npt.assert_array_almost_equal(model.rhs_cst_vector, expected_cst_vector)

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

    def test_build_parametrized(self, name, flux_kwargs, clearance_kwargs, dose, y_init, expected):
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

        
    
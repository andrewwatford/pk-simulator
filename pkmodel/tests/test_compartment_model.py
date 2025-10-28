import numpy.testing as npt
import pytest
import logging
import pkmodel as pk
from pkmodel.builtin_fluxes import constant_dose
import numpy as np

@pytest.fixture()
def cmodel_1():
    """
    Fixture for a simple two-compartment model.
    No fluxes, clearances, or dosages added.
    """
    return pk.CompartmentModel(
        compartment_names=['central', 'peripheral'],
        compartment_volumes=[22, 7]
    )

class TestCompartmentModel:
    """
    Tests the CompartmentModel class.
    """
    def test_create(self, cmodel_1):
        """
        Tests CompartmentModel creation.
        """
        # Check if attributes are stored correctly in a model object
        assert cmodel_1.compartment_names == ["central", "peripheral"]
        assert cmodel_1.compartment_volumes == [22, 7]

    def test_create_with_invalid_inputs(self):
        """
        Tests CompartmentModel creation with mismatched lengths of names and volumes.
        """

        # Test with mismatched lengths of names and volumes
        with pytest.raises(ValueError):
            pk.CompartmentModel(
                compartment_names   = ['central','peripheral'],
                compartment_volumes = [1,2,3])
            
    def test_add_flux_invalid_rate_law(self, cmodel_1):
        """
        Tests that adding a flux with an invalid rate law raises a NotImplementedError.
        """
        
        # Try to add a flux with an invalid rate law
        with pytest.raises(NotImplementedError):
            cmodel_1.add_flux(
                from_compartment = 'central',
                to_compartment   = 'peripheral',
                rate_constant    = 1,
                rate_law         = 'invalid_rate_law')
            
    def test_add_clearance_invalid_rate_law(self, cmodel_1):
        """
        Tests that adding a clearance with an invalid rate law raises a NotImplementedError.
        """     
        
        # Try to add a clearance with an invalid rate law
        with pytest.raises(NotImplementedError):
            cmodel_1.add_clearance(
                from_compartment = 'central',
                rate_constant    = 1,
                rate_law         = 'invalid_rate_law')
            
    def test_add_flux_invalid_nature(self, cmodel_1):
        """
        Tests that adding a flux with an invalid nature raises a NotImplementedError.
        """ 
        
        # Try to add a flux with an invalid nature
        with pytest.raises(NotImplementedError):
            cmodel_1.add_flux(
                from_compartment = 'central',
                to_compartment   = 'peripheral',
                rate_constant    = 1,
                rate_law         = 'first',
                nature           = 'invalid_nature')

    @pytest.mark.parametrize(
        "compartment_names, compartment_volumes, flux_dict_list, clearance_dict_list, expected_matrix, expected_cst_vector",
        [
            # Test case 1: two compartments, one first-order diffusive flux
            (   
                ['central', 'peripheral'],
                [22, 7],
                [{
                    'from_compartment': 'central',
                    'to_compartment': 'peripheral',
                    'rate_constant': 1,
                    'rate_law': 'first',
                    'nature': 'diffusive'}],
                [],
                [
                [-1/22, 1/7],
                [1/22, -1/7]
                ],
                [0, 0]
            ),
            # Test case 2: two compartments, one first-order clearance
            (   
                ['central', 'peripheral'],
                [22, 7],
                [],
                [{
                    'from_compartment': 'central',
                    'rate_constant': 2,
                    'rate_law': 'first'
                }],
                [
                [-2/22, 0],
                [0, 0]
                ],
                [0, 0]
            ),
            # Test case 3: two compartments, one first-order one-way flux
            (   
                ['central', 'peripheral'],
                [22, 7],
                [{
                    'from_compartment': 'central',
                    'to_compartment': 'peripheral',
                    'rate_constant': 1,
                    'rate_law': 'first',
                    'nature': 'one-way'}],
                [],
                [
                [-1/22, 0],
                [1/22, 0]
                ],
                [0, 0]
            ),
            # Test case 4: two compartments, one zero-order flux
            (   
                ['central', 'peripheral'],
                [22, 7],
                [{
                    'from_compartment': 'central',
                    'to_compartment': 'peripheral',
                    'rate_constant': 1,
                    'rate_law': 'zero'}],
                [],
                [
                [0, 0],
                [0, 0]
                ],
                [-1, 1]
            ),
            # Test case 5: two compartments, one zero-order clearance
            (   
                ['central', 'peripheral'],
                [22, 7],
                [],
                [{
                    'from_compartment': 'central',
                    'rate_constant': 2,
                    'rate_law': 'zero'
                }],
                [
                [0, 0],
                [0, 0]
                ],
                [-2, 0]
            ),
            # Test case 6: intravenous bolus model in instructions (dose = 0)
            # Has one first-order diffusive flux between central and peripheral compartments
            # And one first-order clearance from central compartment
            (   
                ['central', 'peripheral'],
                [22, 7],
                [{
                    'from_compartment': 'central',
                    'to_compartment': 'peripheral',
                    'rate_constant': 3,
                    'rate_law': 'first',
                    'nature': 'diffusive'
                }],
                [{
                    'from_compartment': 'central',
                    'rate_constant': 5,
                    'rate_law': 'first'
                }],
                [
                [-8/22, 3/7],
                [3/22, -3/7]
                ],
                [0, 0]
            ),
            # Test case 7: subcutaneous dosing model in instructions. Dosage = 0
            (   
                ['absorber','central', 'peripheral'],
                [5, 22, 7],
                [{
                    'from_compartment': 'absorber',
                    'to_compartment': 'central',
                    'rate_constant': 1,
                    'rate_law': 'first',
                    'nature': 'one-way'
                }, 
                {
                    'from_compartment': 'central',
                    'to_compartment': 'peripheral',
                    'rate_constant': 5,
                    'rate_law': 'first',
                    'nature': 'diffusive'
                }],
                [{
                    'from_compartment': 'central',
                    'rate_constant': 3,
                    'rate_law': 'first'
                }],
                [
                [-1/5, 0, 0],
                [1/5, -8/22, 5/7], 
                [0, 5/22, -5/7]
                ],
                [0, 0, 0]
            ),
            # Test case 8: no fluxes or clearances
            (   
                ['absorber','central', 'peripheral'],
                [5, 22, 7],
                [],
                [],
                [
                [0, 0, 0],
                [0, 0, 0], 
                [0, 0, 0]
                ],
                [0, 0, 0]
            ),
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

        
    
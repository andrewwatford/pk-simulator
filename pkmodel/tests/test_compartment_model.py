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
        "comp_dict, flux_dict, clearance_dict, expected_matrix, expected_cst_vector",
        [
            # Test case 1: two compartments, one first-order diffusive flux
            (   
                {
                    "central":    22.0,
                    "peripheral": 7.0
                },

                {
                    "test_flux": {
                    "source":"central",
                    "dest": "peripheral",
                    "rate_constant": 1.0,
                    "nature":"bidirectional",
                    "rate_law":"first"
                    }
                },
                None,
                [
                    [-1/22, 1/7],
                    [1/22, -1/7]
                ],
                [0, 0]
            ),
            # Test case 2: two compartments, one first-order clearance
            (   
                {
                    "central":    22.0,
                    "peripheral": 7.0
                },
                None,
                {
                    "central_clearance":{
                        "source":"central",
                        "rate_constant": 2.0,
                        "rate_law":"first"
                    }
                },
                [
                    [-2/22, 0],
                    [0, 0]
                ],
                [0, 0]
            ),
            # Test case 3: two compartments, one first-order one-way flux
            (   
                {
                    "central":    22.0,
                    "peripheral": 7.0
                },
                {
                    "test_flux": {
                    "source":"central",
                    "dest": "peripheral",
                    "rate_constant": 1.0,
                    "nature":"unidirectional",
                    "rate_law":"first"
                    }
                },
                None,
                [
                    [-1/22, 0],
                    [1/22,  0]
                ],
                [0, 0]
            ),
            # Test case 4: two compartments, one zero-order flux
            (   
                {
                    "central":    22.0,
                    "peripheral": 7.0
                },
                {
                    "test_flux": {
                    "source":"central",
                    "dest": "peripheral",
                    "rate_constant": 1.0,
                    "nature":"unidirectional",
                    "rate_law":"zero"
                    }
                },
                None,
                [
                    [0, 0],
                    [0, 0]
                ],
                [-1, 1]
            ),
            # Test case 5: two compartments, one zero-order clearance
            (   
                {
                    "central":    22.0,
                    "peripheral": 7.0
                },
                None,
                {
                    "central_clearance":{
                        "source":"central",
                        "rate_constant": 2.0,
                        "rate_law":"zero"
                    }
                },
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
                {   
                    "central":    22.0,
                    "peripheral": 7.0
                },
                {
                    "test_flux": {
                    "source":"central",
                    "dest": "peripheral",
                    "rate_constant": 3.0,
                    "nature":"bidirectional",
                    "rate_law":"first"
                    }
                },
                {
                    "central_clearance":{
                        "source":"central",
                        "rate_constant": 5.0,
                        "rate_law":"first"
                    }
                },
                [
                    [-8/22, 3/7],
                    [3/22, -3/7]
                ],
                [0, 0]
            ),
            # Test case 7: subcutaneous dosing model in instructions. Dosage = 0
            (   
                {
                    "absorber": 5.0,
                    "central": 22.0,
                    "peripheral": 7.0
                },

                {
                    "c_p": {
                        "source":"central",
                        "dest": "peripheral",
                        "rate_constant": 5.0,
                        "nature":"bidirectional",
                        "rate_law":"first"
                    },
                    "a_c": {
                        "source":"absorber",
                        "dest": "central",
                        "rate_constant": 1.0,
                        "nature":"unidirectional",
                        "rate_law":"first"
                    }    
                },

                {
                    "central_clearance":{
                        "source":"central",
                        "rate_constant": 3.0,
                        "rate_law":"first"
                    }
                },
                [
                    [-1/5, 0, 0],
                    [1/5, -8/22, 5/7], 
                    [0, 5/22, -5/7]
                ],
                [0, 0, 0]
            ),
            # Test case 8: no fluxes or clearances
            (   
                {
                    "absorber": 5.0,
                    "central": 22.0,
                    "peripheral": 7.0
            
                },
                None,
                None,
                [
                    [0, 0, 0],
                    [0, 0, 0], 
                    [0, 0, 0]
                ],
                [0, 0, 0]
            ),
        ]
    )
    def test_matrix_and_constant_vector(self, comp_dict, flux_dict, clearance_dict, expected_matrix, expected_cst_vector):
        """
        Tests that the RHS matrix and constant vector are correctly constructed
        after various combinations of zero- or first-order fluxes and clearances are added
        """      

        config = {
            "compartments": comp_dict,
            "fluxes": flux_dict,
            "clearances": clearance_dict,
            "dosages": None
        }

        model = pk.CompartmentModel.from_config(config)
        model.build_linear_rhs()

        npt.assert_array_almost_equal(model.A, expected_matrix)
        npt.assert_array_almost_equal(model.b, expected_cst_vector)

    @pytest.mark.parametrize(
    "name, flux_dict, clearance_dict, dose, y_init, expected" ,
        [
            (
                    "first_order_flux_clearance_dose",
                {
                    "test_flux": {
                    "source":"central",
                    "dest": "peripheral",
                    "rate_constant": 5.0,
                    "nature":"bidirectional",
                    "rate_law":"first"
                    }
                },
                {
                    "central_clearance":{
                        "source":"central",
                        "rate_constant": 5.0,
                        "rate_law":"first"
                    }
                },
                1,  # Dose into the central compartment                    
                np.array([22.0, 7.0]),
                np.array([-4.0, 0.0])
            ),

            (
                    "first_order_flux_dose",
                {
                    "test_flux": {
                    "source":"central",
                    "dest": "peripheral",
                    "rate_constant": 5.0,
                    "nature":"bidirectional",
                    "rate_law":"first"
                    }
                },
                None,
                1,  # Dose into the central compartment                    
                np.array([22.0, 7.0]),
                np.array([1.0, 0.0])
            ),
        ]
    )

    def test_build_parametrized(self, name, flux_dict, clearance_dict, dose, y_init, expected):
        logging.info(f"Testing the model {name}")
        config = {
            "compartments": {
                "central":    22.0,
                "peripheral": 7.0,
            },
            "fluxes": flux_dict,
            "clearances": clearance_dict,
            "dosages": {
                "central_dosage":{
                    "dest":"central",
                    "regime":"constant",
                    "rate_constant": dose,
                }
            }
        }

        model = pk.CompartmentModel.from_config(config)

        # Build the model
        model.build_linear_rhs()

        assert model.model_built == True
        assert model.model_changed_since_last_build == False

        # Testing if the RHS works with some sample numbers
        result = model.rhs(5, np.array(y_init))
        assert result == pytest.approx(np.asarray(expected), abs=1e-8)

        
    
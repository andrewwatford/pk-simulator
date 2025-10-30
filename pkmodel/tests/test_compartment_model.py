import numpy.testing as npt
import pytest
import logging
import pkmodel as pk
from pkmodel.builtin_fluxes import constant_dose
import numpy as np

@pytest.fixture()
def config_1():
    """
    Fixture for a config file.
    """
    config = {

        "compartments": {
            "central":    22.0,
            "peripheral": 7.0,
        },

        "fluxes": {
            "c_p": {
                "source":"central",
                "dest": "peripheral",
                "rate_constant": 5.0,
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
                "regime":"constant",
                "rate_constant": 1.0,
            }
        }
    }

    return config


@pytest.fixture()
def cmodel_1(config_1):
    """
    Fixture for a CompartmentModel instance.
    """
    cmodel = pk.CompartmentModel.from_config(config_1)
    return cmodel


class TestCompartmentModel:
    """
    Tests the CompartmentModel class.
    """
    def test_create(self, cmodel_1):
        """
        Tests CompartmentModel creation.
        """
        from collections import OrderedDict

        assert isinstance(cmodel_1, pk.CompartmentModel)
        assert cmodel_1.model_built == False
        assert cmodel_1.model_changed_since_last_build == True
        assert isinstance(cmodel_1.compartments, OrderedDict)
        assert isinstance(cmodel_1.fluxes, OrderedDict)
        assert isinstance(cmodel_1.clearances, OrderedDict)
        assert isinstance(cmodel_1.dosages, OrderedDict)
        assert len(cmodel_1.compartments) == 2
        assert len(cmodel_1.fluxes) == 1
        assert len(cmodel_1.clearances) == 1
        assert len(cmodel_1.dosages) == 1
        # etc. Could do more

    def test_add_compartment(self, cmodel_1):
        """ 
        Check .add_compartment method works correctly
        """
        cmodel_1.model_changed_since_last_build = False
        cmodel_1.add_compartment(pk.Compartment('peripheral 2', 30))
        assert 'peripheral 2' in cmodel_1.compartments.keys()
        assert cmodel_1.compartments['peripheral 2'].volume == 30
        assert cmodel_1.model_changed_since_last_build == True

    @pytest.mark.parametrize(
            "rate_law_to_use, nature_to_use",
            [('first', 'bidirectional'), ('first', 'unidirectional'), ('zero', 'bidirectional')]
    )
    def test_add_flux(self, cmodel_1, rate_law_to_use, nature_to_use):
        """
        Check .add_flux method works correctly
        """
        cmodel_1.model_changed_since_last_build = False
        cmodel_1.add_flux(pk.Flux(id='test_flux',
                            source = cmodel_1.compartments['central'],
                            dest = cmodel_1.compartments['peripheral'],
                            rate_constant = 1,
                            rate_law = rate_law_to_use,
                            nature = nature_to_use))
        assert 'test_flux' in cmodel_1.fluxes.keys()
        assert cmodel_1.fluxes['test_flux'].rate_law == rate_law_to_use
        assert cmodel_1.fluxes['test_flux'].nature == nature_to_use
        assert cmodel_1.fluxes['test_flux'].source.id == 'central'
        assert cmodel_1.fluxes['test_flux'].dest.id == 'peripheral'
        assert cmodel_1.fluxes['test_flux'].rate_constant == 1
        assert cmodel_1.model_changed_since_last_build == True

    @pytest.mark.parametrize(
            "rate_law_to_use",
            ['first', 'zero']
    )
    def test_add_clearance(self, cmodel_1, rate_law_to_use):
        """
        Check .add_clearance method works correctly
        """
        cmodel_1.model_changed_since_last_build = False
        cmodel_1.add_clearance(pk.Clearance(
                id='test_clearance',
                source = cmodel_1.compartments['central'],
                rate_constant = 3,
                rate_law = rate_law_to_use
            ))
        assert 'test_clearance' in cmodel_1.clearances.keys()
        assert cmodel_1.clearances['test_clearance'].source == cmodel_1.compartments['central']
        assert cmodel_1.clearances['test_clearance'].rate_constant == 3
        assert cmodel_1.clearances['test_clearance'].rate_law == rate_law_to_use
        assert cmodel_1.model_changed_since_last_build == True

    @pytest.mark.parametrize(
            "regime_to_use, rate_constant_to_use, dosage_func_to_use",
            [('constant',0, None), 
             ('constant', 10, None),
             ('custom',0,lambda x: x**2)
             ]
    )
    def test_add_dosage(self, cmodel_1, regime_to_use, rate_constant_to_use, dosage_func_to_use):
        """
        Check .add_dosage method works correctly
        """
        cmodel_1.model_changed_since_last_build = False
        cmodel_1.add_dosage(pk.Dosage(
            id='test_dosage',
            dest=cmodel_1.compartments['central'],
            regime=regime_to_use,
            rate_constant=rate_constant_to_use,
            dosage_func=dosage_func_to_use
        ))
        assert 'test_dosage' in cmodel_1.dosages.keys()
        assert cmodel_1.dosages['test_dosage'].dest == cmodel_1.compartments['central']
        assert cmodel_1.dosages['test_dosage'].regime==regime_to_use
        assert cmodel_1.dosages['test_dosage'].rate_constant==rate_constant_to_use
        assert cmodel_1.dosages['test_dosage'].dosage_func==dosage_func_to_use
        assert cmodel_1.model_changed_since_last_build == True

    # To do:
    # error when you add an existing compartment
    # error when you add an existing flux
    # error when you add a flux to a department that doesn't exist
    # same as above but with clearances
    # same as above but with dosages

    def test_add_flux_invalid_rate_law(self, cmodel_1):
        """
        Tests that adding a flux with an invalid rate law raises a ValueError.
        """
        # Try to add a flux with an invalid rate law
        with pytest.raises(ValueError):
            cmodel_1.add_flux(pk.Flux(id='test',
                            source = cmodel_1.compartments['central'],
                            dest = cmodel_1.compartments['peripheral'],
                            rate_constant = 1,
                            rate_law = "invalid",
                            nature = "unidirectional"))
                    
    def test_add_clearance_invalid_rate_law(self, cmodel_1):
        """
        Tests that adding a clearance with an invalid rate law raises a ValueError.
        """     
        
        # Try to add a clearance with an invalid rate law
        with pytest.raises(ValueError):
            cmodel_1.add_clearance(pk.Clearance(
                id='test',
                source = cmodel_1.compartments['central'],
                rate_constant = 3,
                rate_law = 'invalid'
            ))
            
    def test_add_flux_invalid_nature(self, cmodel_1):
        """
        Tests that adding a flux with an invalid nature raises a ValueError.
        """ 
        
        # Try to add a flux with an invalid nature
        with pytest.raises(ValueError):
            cmodel_1.add_flux(pk.Flux(id='test',
                            source = cmodel_1.compartments['central'],
                            dest = cmodel_1.compartments['peripheral'],
                            rate_constant = 1,
                            rate_law = "invalid",
                            nature = "unidirectional"))

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

        
    
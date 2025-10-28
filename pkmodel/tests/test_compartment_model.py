import numpy.testing as npt
import pytest
import pkmodel as pk

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
        npt.assert_array_equal(model.compartment_names,
            ["central", "peripheral"])
        npt.assert_array_equal(model.compartment_volumes,
            [22, 7])
        
    def test_create_with_invalid_inputs(self):
        """
        Tests CompartmentModel creation with mismatched lengths of names and volumes.
        """

        # Test with mismatched lengths of names and volumes
        with pytest.raises(ValueError):
            pk.CompartmentModel(
                compartment_names   = ['central','peripheral'],
                compartment_volumes = [1,2,3])
    
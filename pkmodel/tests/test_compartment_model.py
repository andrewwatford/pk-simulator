import pytest
import pkmodel as pk

@pytest.fixture()
def test_model():
    pass


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
        
    
import numpy.testing as npt
import pkmodel as pk
from pkmodel.builtin_fluxes import constant_dose
import numpy as np
import pytest

@pytest.mark.parametrize(
    "r_D, r_C, ic",
    [
        (0.1, 1, 0),
        (10, 0.5, 0),
        (5, 0.2, 10),
        (20, 1.0, 5),
    ])
class TestExamples:
    """
    Tests simple examples using the library.
    """
    def test_1d_example(self, r_D, r_C, ic):
        """
        Tests a model with one compartment, constant dosing and first-order clearance.
        """
        # Create a one-compartment model
        model = pk.CompartmentModel(
            compartment_names   = ['central'],
            compartment_volumes = [1])  # Volume
        
        # Add a constant dose flux into the central compartment
        model.add_dosage(
            compartment_name   = 'central',
            dosage_func        = constant_dose(r_D))
        
        # Add a first-order elimination clearance from the central compartment
        model.add_clearance(
            from_compartment = 'central',
            rate_constant    = r_C,  # Clearance rate
            rate_law         = 'first')
        
        # Simulate the model
        t_span = [0, 10]
        y0 = [ic]
        results = model.run(t_span, y0)
        time_points = results.time.data

        # Define the expected result
        expected = ic * np.exp(- r_C * time_points) + (r_D / r_C) * (1 - np.exp(- r_C * time_points))
        
        # Extract mass in the central compartment as array
        central_mass = results['central'].data
        
        # Check that the simulated results match the expected results, to 2 decimal places (any more and numerical errors creep in)
        npt.assert_array_almost_equal(central_mass, expected, decimal = 2)
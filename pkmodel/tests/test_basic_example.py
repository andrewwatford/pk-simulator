import numpy.testing as npt
import pkmodel as pk
from pkmodel.builtin_fluxes import constant_dose
import numpy as np

class TestExamples:
    """
    Tests simple examples using the library.
    """
    def test_1d_example(self):
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
            dosage_func        = constant_dose(1))
        
        # Add a first-order elimination clearance from the central compartment
        model.add_clearance(
            from_compartment = 'central',
            rate_constant    = 0.1,  # Clearance rate
            rate_law         = 'first')
        
        # Simulate the model
        t_span = [0, 48]
        y0 = [0]
        time_points = np.linspace(*t_span, 1000)  # Time points
        results = model.run(t_span, y0, time_points)

        # Define the expected result
        expected = 10 * (1 - np.exp(-0.1 * time_points))
        
        # Extract mass in the central compartment as array
        central_mass = results['central'].data
        
        # Check that the simulated results match the expected results, to 3 decimal places (any more and numerical errors creep in)
        npt.assert_array_almost_equal(central_mass, expected, decimal = 3)
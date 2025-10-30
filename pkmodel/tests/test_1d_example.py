import numpy.testing as npt
import pkmodel as pk
from pkmodel.builtin_fluxes import constant_dose
import numpy as np
from scipy.integrate import quad
import pytest

@pytest.mark.parametrize(
    "r_C, ic, dosage_func",
    [
        (1, 0, constant_dose(0.1)),
        (0.5, 0, constant_dose(10)),
        (0.2, 10, constant_dose(5)),
        (1.0, 5, constant_dose(0)),
        (1.0, 5, lambda t: 10 * np.cos(np.pi * t / 3)**10)
    ])
class Test1dExamples:
    """
    Tests simple examples using the library.
    """
    def test_1d_example(self, r_C, ic, dosage_func):
        """
        Tests a model with one compartment and first-order clearance.
        """
        # Create the compartment, clearance, and dosage
        central = pk.Compartment(id='central', volume=1)
        clearance = pk.Clearance(
            id='cl',
            source=central,
            rate_constant=r_C,
            rate_law='first')
        dosage = pk.Dosage(
            id='dose',
            dest=central,
            regime='custom',
            dosage_func=dosage_func)
        # Create the model and add components
        model = pk.CompartmentModel()
        model.add_compartment(central)
        model.add_clearance(clearance)
        model.add_dosage(dosage)
        # Built
        model.build_linear_rhs()
        # Simulate the model
        t_span = [0, 10]
        y0 = [ic]
        results = model.run(t_span, y0)
        time_points = results.time.data

        # Define the expected result
        def sol_func(t):
            integral, _ = quad(lambda s: np.exp(r_C * s) * dosage_func(s), 0, t)
            return np.exp(- r_C * t) * (ic + integral)
        expected = [sol_func(t) for t in time_points]
        
        # Extract mass in the central compartment as array
        central_mass = results['central'].data
        
        # Check that the simulated results match the expected results, to 2 decimal places (any more and numerical errors creep in)
        npt.assert_allclose(central_mass, expected, rtol = 1e-2)

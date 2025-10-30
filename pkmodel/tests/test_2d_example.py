import numpy.testing as npt
import pkmodel as pk
from pkmodel.builtin_fluxes import constant_dose
import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad_vec
import pytest

@pytest.mark.parametrize(
    "k, r1, r2, ic, d1, d2",
    [
        (1, 1, 1, np.array([0, 0]), constant_dose(0.1), constant_dose(0.2)),
        (1, 1, 1, np.array([0, 0]), constant_dose(0.1), lambda t: 10 * np.cos(np.pi * t / 3)**10),
        (1, 0, 1, np.array([1, 0]), constant_dose(0.1), lambda t: 10 * np.cos(np.pi * t / 3)**10),
        (0, 0, 1, np.array([1, 0]), constant_dose(0.1), lambda t: 10 * np.cos(np.pi * t / 3)**10)
    ])
class Test2dExamples:
    """
    Tests simple examples using the library.
    """
    def test_2d_example(self, k, r1, r2, ic, d1, d2):
        """
        Tests a model with two compartments, each with first-order clearance.
        """
        # Create a two-compartment model
        model = pk.CompartmentModel(
            compartment_names   = ['c1', 'c2'],
            compartment_volumes = [1, 1])  # Volume
        
        # Add the flux between compartments
        model.add_flux(
            from_compartment = 'c1',
            to_compartment   = 'c2',
            rate_constant    = k,
            rate_law         = 'first'
        )

        # Add the two clearances
        model.add_clearance(
            from_compartment = 'c1',
            rate_constant    = r1,
            rate_law         = 'first'
        )
        model.add_clearance(
            from_compartment = 'c2',
            rate_constant    = r2,
            rate_law         = 'first'
        )

        # Add the two dosages
        model.add_dosage(
            compartment_name   = 'c1',
            dosage_func        = d1)
        model.add_dosage(
            compartment_name   = 'c2',
            dosage_func        = d2)
        
        # Simulate the model
        t_span = [0, 10]
        results = model.run(t_span, ic)
        time_points = results.time.data

        # Define the expected result
        def sol_func(t):
            matrix = np.array([[-(k + r1), k],
                               [k, -(k + r2)]])
            mexp_func = lambda s: expm(matrix * s)
            dosage_func = lambda s: np.array([d1(s), d2(s)])
            integral, _ = quad_vec(lambda s: mexp_func(-s) @ dosage_func(s), 0, t)
            return mexp_func(t) @ (ic + integral)
        expected = np.array([sol_func(t) for t in time_points])
        
        # Extract mass in the two compartment as array
        c1_mass = results['c1'].data
        c2_mass = results['c2'].data
        
        # Check that the simulated results match the expected results, to 2 decimal places (any more and numerical errors creep in)
        npt.assert_allclose(c1_mass, expected[:,0], rtol = 1e-2)
        npt.assert_allclose(c2_mass, expected[:,1], rtol = 1e-2)

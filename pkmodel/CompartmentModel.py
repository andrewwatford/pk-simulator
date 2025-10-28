import warnings
from typing import Sequence
#from Compartment import Compartment
#from Flux import Flux
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


def combine_functions(*funcs):
    """Combine multiple scalar time-dependent functions into one vector-valued function.

    ### Args
        *funcs: Callable that accept a simnle float (time) and return a scalar.

    ### Returns:
        A function f(t) that returns a 1-D array with the result of each function.
    """
    def vector_func(x):
        return np.array([f(x) for f in funcs])
    return vector_func

class CompartmentModel:
    """Represents the organism as multiple compartments, with fluxes within and in an out of the system

    ### Attributes:
        compartment_names (list[str]): Names of compartments.
        compartment_volumes (list[float]): Volume of each compartment.
        num_compartments (int): Number of compartments.
        rhs_matrix (np.ndarray): Coefficient matrix for linear contribution to d/dt x (shape n x n).
        rhs_cst_vector (np.ndarray): Constant vector contribution to RHS (shape n).
        dosage_added (list[bool]): Flags whether a dosage function was added per compartment.
        dosage_lst (list[Callable[[float], float]]): Per-compartment dosage functions.
    """
    def __init__(self, compartment_names:(Sequence[str]), compartment_volumes:(Sequence[float])):
        """Initialise the compartment model.

        ### Args:
            compartment_names: Sequence of compartment names.
            compartment_volumes: Sequence of compartment volumes (must match names length).
        """
        if len(compartment_names) != len(compartment_volumes):
            raise ValueError("compartment_names and compartment_volumes must have the same length")
        
        self.compartment_names = compartment_names
        self.compartment_volumes = compartment_volumes

        self.num_compartments = len(self.compartment_names)
        self.rhs_matrix = np.zeros((self.num_compartments, self.num_compartments)) # RHS coefficient matrix
        self.rhs_cst_vector = np.zeros(self.num_compartments) # Constant vector for RHS 
        self.dosage_added = [False for c in compartment_names] # Keeps track of which compartments have a dosage specified
        self.dosage_lst = [lambda t: 0 for c in compartment_names] # Stores all the dosage functions

    def add_flux(self, from_compartment, to_compartment, rate_constant, rate_law='first'):
        """Add a flux between two compartments

        ### Args:
            from_compartment: Source compartment name.
            to_compartment: Destination compartment name.
            rate_constant: Rate constant (units depend on model).
            rate_law: 'first' for first order (proportional to concentration) or 'zero' for zeroth (constant).
        """
        source_idx = self.compartment_names.index(from_compartment)
        source_volume = self.compartment_volumes[source_idx]
        dest_idx = self.compartment_names.index(to_compartment)
        dest_volume = self.compartment_volumes[dest_idx]

        if rate_law == 'first':
            self.rhs_matrix[source_idx, source_idx] += - rate_constant / source_volume
            self.rhs_matrix[source_idx, dest_idx]   += + rate_constant / dest_volume
            self.rhs_matrix[dest_idx, source_idx]   += + rate_constant / source_volume
            self.rhs_matrix[dest_idx, dest_idx]     += - rate_constant / dest_volume
        elif rate_law == 'zero':
            self.rhs_cst_vector[source_idx] += - rate_constant
            self.rhs_cst_vector[dest_idx]   += + rate_constant
        else:
            raise NotImplementedError("Only first or zeroth order fluxes are supported!")
        
    def add_clearance(self, from_compartment, rate_constant, rate_law='first'):
        """Add clearance from a compartment (out of system).

        Args:
            from_compartment: Compartment name where clearance occurs.
            rate_constant: Rate constant (first-order) or constant (zero-order).
            rate_law: 'first' or 'zero'.
        """

        source_idx = self.compartment_names.index(from_compartment)
        source_volume = self.compartment_volumes[source_idx]

        if rate_law == 'first':
            self.rhs_matrix[source_idx, source_idx] += - rate_constant / source_volume
        elif rate_law == 'zero':
            self.rhs_cst_vector[source_idx] += - rate_constant
        else:
            raise NotImplementedError("Wnly first or zeroth order clearances are supported!")
    
    def add_dosage(self, compartment_name, dosage_func):
        """Add a time-dependent dosage function to a compartment.

        Args:
            compartment_name: Which compartment receives dosage.
            dosage_func: Callable f(t) returning the dosage rate at time t.

        Notes:
            Only one dosage function per compartment is stored; subsequent calls overwrite.
        """
         
        compartment_idx = self.compartment_names.index(compartment_name)

        if self.dosage_added[compartment_idx]:
            warnings.warn("Only one dosage regime per compartment supported. Overwriting with most recent dosage.")

        self.dosage_lst[compartment_idx] = dosage_func
        self.dosage_added[compartment_idx] = True


    def build(self):
        """Build and ODE system to be solved with scipy.integrate.solve_ivp"""
        dosage_func_vector = combine_functions(*self.dosage_lst)
        def rhs(t, y):
            dydt = self.rhs_matrix @ y + self.rhs_cst_vector + dosage_func_vector(t)
            return dydt
        return rhs
    
    def run(self, t_span, y0, t_eval=None):
        rhs = self.build()
        sol = solve_ivp(rhs, t_span, y0, t_eval=t_eval, vectorized=False)
        return sol
    
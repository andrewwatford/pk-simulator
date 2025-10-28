import warnings
from typing import Sequence, Callable
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
import logging


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
        model_built (bool): Flag to track if the model has been built.
        model_changed_since_last_build (bool): Flag to track if changes have been made to the model
            in which case it needs to be rebult.
        rhs (Callable[[float, np.ndarray], np.ndarray]): RHS function that computes dydt = f(t, y).
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

        self.model_built = False 
        self.model_changed_since_last_build = True

    def add_flux(self, from_compartment: str, to_compartment: str, rate_constant: float, rate_law: str = "first"):
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
        
        self.model_changed_since_last_build = True # System state has changed, model needs to be rebuilt before running
        
    def add_clearance(self, from_compartment: str, rate_constant: float, rate_law: str = "first"):
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
            raise NotImplementedError("Only first or zeroth order clearances are supported!")
        
        self.model_changed_since_last_build = True # System state has changed, model needs to be rebuilt before running
    
    def add_dosage(self, compartment_name: str, dosage_func: Callable[[float], float]):
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

        self.model_changed_since_last_build = True # System state has changed, model needs to be rebuilt before running

    def build(self):
        """Build and ODE system to be solved with scipy.integrate.solve_ivp
        
        ### Returns:
        Callable[[float, np.ndarray], np.ndarray]: RHS function that computes dydt = f(t, y).
        """
        dosage_func_vector = combine_functions(*self.dosage_lst)
        def rhs(t, y):
            dydt = self.rhs_matrix @ y + self.rhs_cst_vector + dosage_func_vector(t)
            return dydt
        
        self.model_changed_since_last_build = False
        self.model_built = False
        return rhs
    
    def run(self, t_span:Sequence[float], y0:Sequence[float], t_eval:Sequence[float]=None):
        if not self.model_built:
            logging.info("No build detected, building the model from scratch...")
            self.rhs = self.build()
        else:
            if not self.model_changed_since_last_build:
                logging.info("No changes detected since last build. Using the existing build.")
            else:
                logging.info("Changes to the model detected since last build. Rebuilding the model...")
                self.rhs = self.build()

        sol = solve_ivp(self.rhs, t_span, y0, t_eval=t_eval, vectorized=True)
        return sol



  





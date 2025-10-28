import warnings
from scipy.integrate import solve_ivp
import numpy as np

def combine_functions(*funcs):
    def vector_func(x):
        return np.array([f(x) for f in funcs])
    return vector_func

class CompartmentModel:

    def __init__(self, compartment_names, compartment_volumes):
        self.compartment_names = compartment_names
        self.compartment_volumes = compartment_volumes
        self.num_compartments = len(self.compartment_names)
        self.rhs_matrix = np.zeros((self.num_compartments, self.num_compartments))
        self.rhs_cst_vector = np.zeros(self.num_compartments)
        self.dosage_added = [False for c in compartment_names]
        self.dosage_lst = [lambda t: 0 for c in compartment_names]

    def add_flux(self, from_compartment, to_compartment, rate_constant, rate_law='first'):
        source_idx = self.compartment_names.index(from_compartment)
        source_volume = self.compartment_volumes[source_idx]
        dest_idx = self.compartment_names.index(to_compartment)
        dest_volume = self.compartment_volumes[dest_idx]
        if rate_law == 'first':
            self.rhs_matrix[source_idx, source_idx] += - rate_constant / source_volume
            self.rhs_matrix[source_idx, dest_idx] += rate_constant / dest_volume
            self.rhs_matrix[dest_idx, source_idx] += rate_constant / source_volume
            self.rhs_matrix[dest_idx, dest_idx] += - rate_constant / dest_volume
        if rate_law == 'zero':
            self.rhs_cst_vector[source_idx] += -rate_constant
            self.rhs_cst_vector[dest_idx] += +rate_constant
        else:
            raise NotImplementedError("We only support first or zero order fluxes!")
        
    def add_clearance(self, from_compartment, rate_constant, rate_law='first'):
        source_idx = self.compartment_names.index(from_compartment)
        source_volume = self.compartment_volumes[source_idx]
        if rate_law == 'first':
            self.rhs_matrix[source_idx, source_idx] += - rate_constant / source_volume
        if rate_law == 'zero':
            self.rhs_cst_vector[source_idx] += -rate_constant
        else:
            raise NotImplementedError("We only support first or zero order clearances!")
    
    def add_dosage(self, compartment_name, dosage_func):
        compartment_idx = self.compartment_names.index(compartment_name)
        if self.dosage_added[compartment_idx]:
            warnings.warn("Only one dosage per compartment. Overwriting with most recent dosage.")
        self.dosage_lst[compartment_idx] = dosage_func
        self.dosage_added[compartment_idx] = True




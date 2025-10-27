from collections import OrderedDict
from Compartment import Compartment
from Flux import Flux

class CompartmentModel:

    def __init__(self, compartment_names, compartment_volumes):
        self.compartments = OrderedDict()
        self.fluxes = []
        for name, vol in zip(compartment_names, compartment_volumes):
            self.compartments[name] = Compartment(name, vol)

    def add_dosage(self, compartment_name, dosage_function):
        self.compartments[compartment_name].set_dosage(dosage_function)

    def add_clearance(self, compartment_name, clearance_function):
        self.compartments[compartment_name].set_clearance(clearance_function)

    def add_flux(self, from_compartment, to_compartment, rate_function):
        self.fluxes.append(Flux(from_compartment, to_compartment, rate_function))
            
        


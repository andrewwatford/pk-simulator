"""
pkmodel is a Python library for creating, solving, and visualising pharmacokinetic (PK) models.
It allows the user to represent an organism as a set of interacting compartments and simulate the movement of substances between them over time. 
"""
# Import version info
from .version_info import VERSION_INT, VERSION  # noqa

# Import main classes
from .CompartmentModel import CompartmentModel, Compartment, Clearance, Flux, Dosage    # noqa
from .builtin_fluxes import constant_dose   # noqa

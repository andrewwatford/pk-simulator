"""pkmodel is a Pharmokinetic modelling library.

It contains functionality for creating, solving, and visualising the solution
of Parmokinetic (PK) models

"""
# Import version info
from .version_info import VERSION_INT, VERSION  # noqa

# Import main classes
from .CompartmentModel import CompartmentModel, Compartment, Clearance, Flux, Dosage    # noqa
from .builtin_fluxes import constant_dose   # noqa

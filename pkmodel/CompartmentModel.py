# Data types and classes
from typing import Sequence, Callable
from collections import OrderedDict
from dataclasses import dataclass
from enum import Enum

# Maths + plotting
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import xarray as xr

# Logging
import warnings
import logging

# Module packages
from .builtin_fluxes import constant_dose

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

@dataclass
class Compartment:
    """Represents a compartment/reservoir within the organism"""
    id: str
    volume: float

@dataclass
class Flux:
    """Represents a flux connecting two compartments"""
    id: str
    source: Compartment
    dest:   Compartment
    rate_constant: float
    rate_law: str = "first"
    nature: str = "unidirectional"
    
    def __post_init__(self):
        if self.rate_law not in ['first', 'zero']:
            raise ValueError(f"Rate law '{self.rate_law}' is not supported! Supported rate laws are 'first' and 'zero'.")
        if self.nature not in ['unidirectional', 'bidirectional']:
            raise ValueError(f"Nature '{self.nature}' is not supported! Supported natures are 'unidirectional' and 'bidirectional'.")

@dataclass 
class Clearance:
    """Represents a compound leaving a compartment"""
    id: str
    source: Compartment
    rate_constant: float
    rate_law: str = 'first'

    def __post_init__(self):
        if self.rate_law not in ['first', 'zero']:
            raise ValueError(f"Rate law '{self.rate_law}' is not supported! Supported rate laws are 'first' and 'zero'.")

@dataclass
class Dosage:
    """Represents a compound entering a compartment"""
    id: str
    dest: Compartment
    regime: str = 'constant'
    rate_constant: float = 0  # This rate constant is used for the constant dosing regime
    dosage_func: Callable = None

    def __post_init__(self):
        if self.regime not in ['constant', 'custom']:
            raise ValueError(f"Dosage regime '{self.regime}' is not supported! Supported regimes are 'constant' and 'custom'.")

class CompartmentModel:
    """Represents the organism as multiple compartments, with fluxes within and in an out of the system

    ### Attributes:
        compartments (OrderedDict[str: Compartment])
        fluxes (OrderedDict[str: Flux])
        dosages (OrderedDict[str: Clearance])
        clearances (OrderedDict[str: Dosage])

        model_built (bool): Flag to track if the model has been built.
        model_changed_since_last_build (bool): Flag to track if changes have been made to the model
            in which case it needs to be rebult.

        rhs (Callable[[float, np.ndarray], np.ndarray]): RHS function that computes dydt = f(t, y).
    """
    def __init__(self):
        """Initialise the compartment model."""

        self.compartments: "OrderedDict[str: Compartment]"  = OrderedDict()
        self.fluxes:       "OrderedDict[str: Flux]"         = OrderedDict()
        self.clearances:   "OrderedDict[str: Clearance]"    = OrderedDict()
        self.dosages:      "OrderedDict[str: Dosage]"       = OrderedDict()
        
        # Flags to check if the model has been built
        self.model_built = False 
        self.model_changed_since_last_build = True

    def add_compartment(self, comp:Compartment):
        if comp.id in self.compartments:
            raise KeyError(f"Compartment with id '{comp.id}' alredy exists in the model!")
        
        self.compartments[comp.id] = comp
        self.model_changed_since_last_build = True 

    def add_flux(self, flux:Flux):
        if flux.id in self.fluxes:
            raise KeyError(f"Flux with id '{flux.id}' alredy exists in the model!")
        
        if flux.source.id not in self.compartments or flux.dest.id not in self.compartments:
            raise KeyError(f"Can't add flux to the model: one or both compartments {flux.source.id} or {flux.dest.id} are not present!")
        
        self.fluxes[flux.id] = flux
        self.model_changed_since_last_build = True

    def add_clearance(self, clr:Clearance):
        if clr.id in self.clearances:
            raise KeyError(f"Clerance with id '{clr.id}' alredy exists in the model!")
        
        if clr.source.id not in self.compartments:
            raise KeyError(f"Can't a clearance to the model: compartment {clr.source.id} is not present")
        
        self.clearances[clr.id] = clr
        self.model_changed_since_last_build = True

    def add_dosage(self, dsg:Dosage):
        if dsg.id in self.dosages:
            raise KeyError(f"Dosage with id '{dsg.id}' alredy exists in the model!")
        
        if dsg.dest.id not in self.compartments:
            raise KeyError(f"Can't a soxage to the model: compartment {dsg.dest.id} is not present")
        
        self.dosages[dsg.id] = dsg
        self.model_changed_since_last_build = True

    @classmethod
    def from_config(cls, config): #TODO check the config against a schema, add support for missing arguments
        model = cls()

        for id, vol in config["compartments"].items():
            model.add_compartment(
                Compartment(
                    id=id,
                    volume=vol
                )
            )

        if config.get("fluxes") != None:     
            for id, flux in config["fluxes"].items():
                model.add_flux(
                    Flux(
                        id=id,
                        source=model.compartments[flux["source"]],
                        dest=model.compartments[flux["dest"]],
                        rate_constant = flux["rate_constant"],
                        nature = flux["nature"],
                        rate_law = flux["rate_law"]
                    )
                )

        if config.get("clearances") != None:
            for id, clr in config["clearances"].items():
                model.add_clearance(
                    Clearance(
                        id=id,
                        source=model.compartments[clr["source"]],
                        rate_constant = clr["rate_constant"],
                        rate_law = clr["rate_law"]
                    )
                )

        if config.get("dosages") != None:
            for id, dsg in config["dosages"].items():
                model.add_dosage(
                    Dosage(
                        id=id,
                        dest=model.compartments[dsg["dest"]],
                        regime=dsg["regime"],
                        rate_constant = dsg["rate_constant"],
                    )
                )

        return model
        
    def build_linear_rhs(self):
        
        n = len(self.compartments)
        A = np.zeros((n, n), dtype = float) # RHS coefficient matrix
        b = np.zeros(n, dtype=float) # Constant vector for RHS 

        # Build a numeric index of compartments
        self.comp_index = {comp.id: i for i, comp in enumerate(self.compartments.values())} # TODO: this may be moved to a separate method in the future

        # Fluxes
        for flux in self.fluxes.values():
            src_idx = self.comp_index[flux.source.id]
            dst_idx   = self.comp_index[flux.dest.id]
            if flux.rate_law == 'first':
                A[src_idx, src_idx] += - flux.rate_constant / flux.source.volume
                A[dst_idx, src_idx] += + flux.rate_constant / flux.source.volume

                if flux.nature == 'bidirectional':
                    A[src_idx, dst_idx]   += + flux.rate_constant / flux.dest.volume
                    A[dst_idx, dst_idx]   += - flux.rate_constant / flux.dest.volume

            elif flux.rate_law == 'zero':
                warnings.warn("Zero order fluxes are supported, but implementation is non-physical.")
                b[src_idx] += - flux.rate_constant
                b[dst_idx] += + flux.rate_constant

        # Clearances
        for clr in self.clearances.values():
            src_idx = self.comp_index[clr.source.id]
            if clr.rate_law == 'first':
                A[src_idx, src_idx] += - clr.rate_constant / clr.source.volume
            elif clr.rate_law == 'zero':
                warnings.warn("Zero order clerances are supported, but implementation is non-physical.")
                b[src_idx] += - clr.rate_constant
            else: 
                raise NotImplementedError("Only first or zero order clearances are supported!")

        # Dosages
        dosage_lst = [lambda t: 0 for c in range(n)] # Stores all the dosage functions
        for dsg in self.dosages.values():
            dst_idx = self.comp_index[dsg.dest.id]
            if dsg.regime == 'constant':      
                dosage_lst[dst_idx] = constant_dose(dsg.rate_constant)
            else: # Custom dosage function
                dosage_lst[dst_idx] = dsg.dosage_func

        d = combine_functions(*dosage_lst)

        # Build the final callable
        def rhs(t, y):
            dydt = A @ y + b + d(t)
            return dydt
        
        self.A = A
        self.b = b
        self.rhs = rhs
        self.model_built = True
        self.model_changed_since_last_build = False
    
    def run(self, t_span:Sequence[float], y0:Sequence[float], t_eval:Sequence[float]=None):
        """Solves (and builds, if self.build() has not been called previously) the ODE system.

        ### Args:
            - t_span: Sequence[float]. A two-element sequence with the start and end times.
            - y0: Sequence[float]. The initial condition for the dynamical system. Order of
                variables matches the order of the compartments specified in CompartmentModel
                construction.
            - t_eval: Sequence[float] (default None). The sequence of time points to solve the
                system on.

        ### Returns:
            - ds: XArray.Dataset. An XArray Dataset object containing the labelled output of the
                simulation.
        """
        if not self.model_built:
            logging.info("No build detected, building the model from scratch...")
            self.build()
        else:
            if not self.model_changed_since_last_build:
                logging.info("No changes detected since last build. Using the existing build.")
            else:
                logging.info("Changes to the model detected since last build. Rebuilding the model...")
                self.build()

        sol = solve_ivp(self.rhs, t_span, y0, t_eval=t_eval, vectorized=False)
        da_dct = {}
        for idx, name in enumerate(self.compartments.keys()):
            da_dct[name] = xr.DataArray(data = sol.y[idx, :], coords = {'time': sol.t})
        ds = xr.Dataset(da_dct)
        return ds
    
    def plot_all(self, ds: xr.Dataset):
        """Plot all compartments from the output dataset.

        ### Args:
            - ds: XArray.Dataset. The output dataset from a model run.

        ### Returns:
            - fig: matplotlib.figure.Figure: The figure object containing the plots.
            - axs: np.ndarray: The array of axes objects for each compartment plot.
        """
        num_compartments = len(self.compartments)
        fig, axs = plt.subplots(num_compartments, 1, sharex=True)
        for idx, name in enumerate(self.compartments.keys()):
            ds[name].plot(ax=axs[idx])
            axs[idx].set_ylabel(f'$q_{{{name}}}$')
            if idx < num_compartments - 1:
                axs[idx].set_xlabel(None)
        axs[-1].set_xlabel('$t$')
        axs[0].set_title('Compartment masses over time')
        return fig, axs

    def construct_graph(self):
        """Construct a NetworkX graph representation of the compartment model.

        ### Returns:
            - g: networkx.MultiDiGraph. A directed multigraph representing the compartment model.
        """
        # Create the empty directed multigraph
        g = nx.MultiDiGraph()
        # Add generic IN and OUT nodes
        g.add_node("IN", subset = "dosages")
        g.add_node("OUT", subset = "clearances")
        # Add compartments as nodes, as well as a generic IN and OUT node for each (TODO: improve representation of these)
        for comp_name, comp in self.compartments.items():
            g.add_node(comp_name, subset = "compartments", **comp.__dict__)
        # Add fluxes, clearances, and dosages as edges
        for flux_name, flux in self.fluxes.items():
            g.add_edge(flux.source.id, flux.dest.id, key=flux_name, **flux.__dict__)
            if flux.nature == "bidirectional":
                g.add_edge(flux.dest.id, flux.source.id, key=flux_name+"_rev", **flux.__dict__)
        for clear_name, clear in self.clearances.items():
            g.add_edge(clear.source.id, "OUT", key=clear_name, **clear.__dict__)
        for dose_name, dose in self.dosages.items():
            g.add_edge("IN", dose.dest.id, key=dose_name, **dose.__dict__)
        self.graph = g
        return self.graph

    def draw_basic_graph(
            self,
            pos: dict = None,
            node_shape: str = "s",
            node_size: int = 3000,
            font_size: int = 10,
            node_color: str = "white",
            edge_color: str = "black",
            linewidths: int = 2,
            arrowsize: int = 20,
            ):
        """Create a basic plot of the compartment model graph.

        ### Returns:
            - fig: matplotlib.figure.Figure. The figure object containing the plot.
            - ax: matplotlib.axes.Axes. The axes object for the plot.
        """
        warnings.warn("This feature is still under development. To produce your own drawing, you may use the output of construct_graph along with NetworkX.")
        if not hasattr(self, 'graph'):
            self.construct_graph()
        g = self.graph
        if pos is None:
            # pos = nx.multipartite_layout(g, align = "horizontal") # TODO: need a better layout
            pos = nx.planar_layout(g) # TODO: need a better layout
        visible_nodes = list(self.compartments.keys())
        visible_node_labels = {node: node.capitalize() for node in visible_nodes}

        fig, ax = plt.subplots()
        nx.draw_networkx_nodes(
            g, 
            pos=pos, 
            nodelist=visible_nodes, 
            node_shape=node_shape,
            node_size=node_size,
            node_color=node_color,
            edgecolors=edge_color,
            linewidths=linewidths
        )
        nx.draw_networkx_labels(
            g, 
            pos=pos, 
            labels=visible_node_labels, 
            font_size=font_size
        )
        nx.draw_networkx_edges(
            g, 
            pos=pos, 
            width=linewidths, 
            arrowsize=arrowsize, 
            node_size=node_size
        )

        # Set margins for the axes so that nodes aren't clipped
        ax = plt.gca()
        ax.margins(0.20)
        plt.axis("off")
        return fig, ax

if __name__ == "__main__":
    central = Compartment(
        id="central", # TODO make it generate an id if it's not provided?
        volume=22
    )
    peripheral = Compartment(
        id="peripheral",
        volume=7
    )

    c_p_flux = Flux(
        id="c_p_flux",
        source=central,
        dest=peripheral,
        rate_constant=5,
        nature="bidirectional",
        rate_law="first"
    )

    central_clr = Clearance(
        id="central_clearance",
        source=central,
        rate_constant=5,
        rate_law="first"
    )

    central_dsg = Dosage(
        id="central_dosage",
        dest=central,
        regime="constant",
        rate_constant=1
    )

    model = CompartmentModel()
    
    model.add_compartment(central)
    model.add_compartment(peripheral)

    model.add_flux(c_p_flux)

    model.add_clearance(central_clr)

    model.add_dosage(central_dsg)

    model.build_linear_rhs()

    import pytest
    y_init = np.array([22.0, 7.0])
    expected = np.array([-4.0, 0.0])
    result = model.rhs(5, np.array(y_init))
    assert result == pytest.approx(np.asarray(expected), abs=1e-8)
    print("Success")
    

    config = {

        "compartments": {
            "central":    22.0,
            "peripheral": 7.0,
        },

        "fluxes": {
            "c_p": {
                "source":"central",
                "dest": "peripheral",
                "rate_constant": 5.0,
                "nature":"bidirectional",
                "rate_law":"first"
            }
        },

        "clearances": {
            "central_clearance":{
                "source":"central",
                "rate_constant": 5.0,
                "rate_law":"first"
            }
        },

        "dosages": {
            "central_dosage":{
                "dest":"central",
                "regime":"constant",
                "rate_constant": 1.0,
            }
        }
    }

    model2 = CompartmentModel.from_config(config)
    model2.build_linear_rhs()
    result = model2.rhs(5, np.array(y_init))
    assert result == pytest.approx(np.asarray(expected), abs=1e-8)
    print("Success")

    config3 = {

        "compartments": {
            "central":    22.0,
            "peripheral": 7.0,
        },

        "fluxes": {
            "c_p": {
                "source":"central",
                "dest": "peripheral",
                "rate_constant": 5.0,
                "nature":"bidirectional",
                "rate_law":"first"
            }
        },

        "clearances": None,

        "dosages": {
            "central_dosage":{
                "dest":"central",
                "regime":"constant",
                "rate_constant": 1.0,
            }
        }
    }

    model3 = CompartmentModel.from_config(config3)
    model3.build_linear_rhs()
    result = model3.rhs(5, np.array(y_init))
    expected=np.array([1.0, 0.0])
    assert result == pytest.approx(np.asarray(expected), abs=1e-8)
    print("Success")

    # TODO - nice print statments for all the classes
    # TODO - check if classes get copied or smth
    # TODO - get a config from a file


    # identity check (preferred)
    if c_p_flux.source is central_clr.source:
        print("They are the same object (identity).")
    else:
        print("Different objects (even if equal by ==).")
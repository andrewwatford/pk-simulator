import pytest
import pkmodel as pk


@pytest.fixture()
def config_1():
    """
    Fixture for a config file.
    """
    config = {

        "compartments": {
            "central": 22.0,
            "peripheral": 7.0,
        },

        "fluxes": {
            "c_p": {
                "source": "central",
                "dest": "peripheral",
                "rate_constant": 5.0,
                "nature": "bidirectional",
                "rate_law": "first"
            }
        },

        "clearances": {
            "central_clearance": {
                "source": "central",
                "rate_constant": 5.0,
                "rate_law": "first"
            }
        },

        "dosages": {
            "central_dosage": {
                "dest": "central",
                "regime": "constant",
                "rate_constant": 1.0,
            }
        }
    }

    return config


@pytest.fixture()
def cmodel_1(config_1):
    """
    Fixture for a CompartmentModel instance.
    """
    cmodel = pk.CompartmentModel.from_config(config_1)
    return cmodel


class TestGraphCreation:
    """
    Test graph creation from CompartmentModel.
    """

    def test_graph_creation(self, cmodel_1):
        """
        Test graph creation from CompartmentModel.
        """
        graph = cmodel_1.construct_graph()

        # Check nodes
        expected_nodes = {'central', 'peripheral', 'central_IN', 'central_OUT', 'peripheral_IN', 'peripheral_OUT'}
        assert set(graph.nodes) == expected_nodes

        # Check edges
        expected_edges = {
            ('central', 'peripheral', 'c_p'),
            ('central', 'central_OUT', 'central_clearance'),
            ('central_IN', 'central', 'central_dosage')
        }
        assert set(graph.edges) == expected_edges

import unittest
import pkmodel as pk


class ComparmentTest(unittest.TestCase):
    """
    Tests the Compartment class.
    """
    def test_create(self):
        """
        Tests Compartment creation.
        """
        compartment = pk.Compartment(
            name = "central",
            volume = 25,
            dosage = None,
            clearance = None
        )
        
        self.assertEqual(compartment.name, "central")

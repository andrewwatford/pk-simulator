class Compartment:

    def __init__(self, name, volume):
        self.name = name
        self.volume = volume
        self.dosage = None
        self.clearance = None

    def set_dosage(self, dosage_function):
        self.dosage = dosage_function

    def set_clearance(self, clearance_function):
        self.clearance = clearance_function
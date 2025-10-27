class Flux:
    def __init__(self, source_compartment_name, dest_compartment_name, rate_function):
        self.source_compartment_name = source_compartment_name
        self.dest_compartment_name = dest_compartment_name
        self.rate_function = rate_function
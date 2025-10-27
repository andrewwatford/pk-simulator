def first_order_flux(c1, c2, k = 0.5):
    """Calculates the first-order kinetics flux law.
    
    Parameters
    ----------
    c1 : float
        Concentration in the source compartment. [mass/volume]
    c2 : float
        Concentration in the destination compartment. [mass/volume]
    k : float
        Rate constant. [volume/time]

    Convention is that flux is positive when going from compartment 1 to compartment 2.
    """
    return k * (c1 - c2)

def constant_dose(t, rate = 0.0):
    """
    Returns a constant dose rate.
    
    Parameters
    ----------
    t : float
        Time. [time]
    rate : float
        Dose rate. [mass/time]
    """
    return rate

def first_order_clearance(c, k = 0.5):
    """
    Calculates the first-order clearance rate.
    
    Parameters
    ----------
    c : float
        Concentration in the compartment. [mass/volume]
    k : float
        Clearance rate constant. [volume/time]
    """
    return k * c
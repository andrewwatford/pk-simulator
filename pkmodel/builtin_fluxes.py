def constant_dose(rate):
    """
    Returns a function that doses at a constant rate.
    
    Parameters
    ----------
    rate : float
        Dose rate. [mass/time]

    Returns
    -------
    dose_func : function
        Function that returns the dose rate at time t. [mass/time]
    """
    def dose_func(t):
        """
        Function that returns the dose rate at time t.

        Parameters
        ----------
        t : float
            Time. [time]

        Returns
        -------
        rate : float
            Dose rate. [mass/time]
        """
        return rate
    return dose_func

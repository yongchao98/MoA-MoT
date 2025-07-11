def get_function_name():
    """
    This function provides the common name for the space-time, double Fourier
    transform of the generalized pair correlation function in nuclear criticality.
    """
    # The generalized pair correlation function describes the correlation of
    # neutron events in space and time.
    # Its double Fourier transform (in space and time) is a key quantity
    # in neutron noise analysis, relating to the power spectrum of the
    # neutron population fluctuations.
    function_name = "Spectral Density Function"
    return function_name

# Print the answer to the user's question.
name = get_function_name()
print(f"In the nuclear criticality community, the space-time, double Fourier transform of the generalized pair correlation function is commonly called the: {name}")

def print_transform_name():
    """
    This script identifies and prints the common name for the space-time,
    double Fourier transform of the generalized pair correlation function
    as used in the nuclear criticality and reactor physics community.
    """
    # The generalized pair correlation function measures the correlation
    # of neutron density fluctuations in space and time.
    # C(r, t) = <delta_n(r_0, t_0) * delta_n(r_0 + r, t_0 + t)>

    # A space-time, double Fourier transform converts this function from the
    # space-time domain (r, t) to the wavevector-frequency domain (k, omega).
    # S(k, omega) = FT[C(r, t)]

    # The resulting function, S(k, omega), represents the fluctuation power
    # spectrum.
    function_name = "Power Spectral Density (PSD)"

    print("In the nuclear criticality community, the space-time, double Fourier transform of the generalized pair correlation function is commonly called the:")
    print(function_name)

if __name__ == "__main__":
    print_transform_name()
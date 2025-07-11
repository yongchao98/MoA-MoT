def print_energy_spectrum():
    """
    This function prints the final expression for the energy spectrum
    for a harmonic oscillator with a quartic perturbation, as predicted
    by the given self-energy diagram.
    """

    # The calculation steps are detailed in the text above.
    # This script formats and prints the final symbolic result.

    print("The self-energy diagram leads to a renormalization of the oscillator frequency.")
    print("The new frequency omega' is approximately omega_0 + (u * hbar) / (8 * m^2 * omega_0^2).")
    print("The predicted energy spectrum E_n, up to an overall constant, is that of a harmonic oscillator with this new frequency:")
    print("E_n = hbar * omega' * (n + 1/2)")
    print("\nSubstituting the expression for omega', the final equation for the energy spectrum is:")
    print("E_n = hbar * (omega_0 + (u * hbar) / (8 * m**2 * omega_0**2)) * (n + 1/2)")
    print("\nExpanding this, we get the final result, with each number shown explicitly:")
    print("E_n = hbar * omega_0 * (n + 1/2) + (u * hbar**2) / (8 * m**2 * omega_0**2) * (n + 1/2)")

print_energy_spectrum()
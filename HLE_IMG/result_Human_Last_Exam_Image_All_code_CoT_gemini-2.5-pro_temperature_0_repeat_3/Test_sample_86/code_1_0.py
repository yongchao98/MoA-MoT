def print_energy_spectrum_equation():
    """
    This function prints the final equation for the energy spectrum,
    outputting each number explicitly as requested.
    The equation is: E_n = E_0 + n * hbar * sqrt(omega_0**2 + (u * hbar) / (4 * m**2 * omega_0))
    """
    print("E_n = E_0 + n * hbar * sqrt(omega_0**", end="")
    print(2, end="")
    print(" + (u * hbar) / (", end="")
    print(4, end="")
    print(" * m**", end="")
    print(2, end="")
    print(" * omega_0))")

print_energy_spectrum_equation()
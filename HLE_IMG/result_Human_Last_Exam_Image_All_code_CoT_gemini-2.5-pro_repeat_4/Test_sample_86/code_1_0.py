def solve_energy_spectrum():
    """
    This function prints the derived formula for the energy spectrum of the
    harmonic oscillator with a quartic perturbation, as predicted by the
    first-order self-energy diagram.
    """
    
    # The coefficients in the final equation as requested.
    coeff_n = 1
    coeff_half = 2
    coeff_denominator = 4
    
    # Constructing the string for the final equation.
    # The equation is E_n = (n + 1/2) * hbar * sqrt(omega_0^2 + (u * hbar) / (4 * m^2 * omega_0))
    equation = (
        "E_n = (n + " + str(coeff_n) + "/" + str(coeff_half) + ") * hbar * "
        "sqrt(omega_0^2 + (u * hbar) / (" + str(coeff_denominator) + " * m^2 * omega_0))"
    )
    
    print("The energy spectrum E_n is given by the following equation:")
    print(equation)

solve_energy_spectrum()
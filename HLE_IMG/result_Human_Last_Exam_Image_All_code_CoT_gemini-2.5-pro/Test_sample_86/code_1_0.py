def solve_harmonic_oscillator_perturbation():
    """
    This function prints the final equation for the energy spectrum of a
    harmonic oscillator with a quartic perturbation, as predicted by the
    self-energy tadpole diagram.
    """

    # The formula for the energy spectrum E_n is derived from the renormalized
    # frequency omega'. The original energy is E_n_0 = (n + 1/2)*hbar*omega_0.
    # The new frequency is omega' = omega_0 + (u*hbar)/(8*m**2*omega_0**2).
    # The new energy spectrum is E_n = (n + 1/2)*hbar*omega'.

    # We construct the string representing the final equation, including all numbers.
    # The numbers in the equation are 1, 2, and 8.
    equation_part1 = "hbar * omega_0 * (n + 1/2)"
    equation_part2 = "(u * hbar**2) / (8 * m**2 * omega_0**2) * (n + 1/2)"

    final_equation = f"E_n = {equation_part1} + {equation_part2}"

    print("The predicted energy spectrum is:")
    print(final_equation)

solve_harmonic_oscillator_perturbation()
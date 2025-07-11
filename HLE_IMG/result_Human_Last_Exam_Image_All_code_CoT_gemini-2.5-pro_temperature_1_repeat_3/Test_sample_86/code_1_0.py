def solve_energy_spectrum():
    """
    This function prints the derived formula for the energy spectrum of a harmonic
    oscillator with a quartic perturbation, based on the tadpole self-energy diagram.
    The formula is constructed and printed piece by piece to show the numbers involved.
    """

    # The formula for the energy levels of a harmonic oscillator is E_n = (n + 1/2) * hbar * w'
    # where w' is the new, renormalized frequency.
    
    # Define the numerical constants in the formula
    n_plus_factor_numerator = 1
    n_plus_factor_denominator = 2
    u_hbar_denominator_coeff = 4
    m_exponent = 2

    # Print the final equation using f-string formatting to highlight the numbers.
    print("The self-energy diagram predicts an energy spectrum corresponding to a renormalized harmonic oscillator.")
    print("The predicted energy spectrum E_n is given by the formula:")
    print(f"E_n = (n + {n_plus_factor_numerator}/{n_plus_factor_denominator}) * hbar * sqrt(w0**2 + (u*hbar) / ({u_hbar_denominator_coeff}*m**{m_exponent}*w0))")

solve_energy_spectrum()
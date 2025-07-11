def solve_energy_spectrum():
    """
    This function calculates and prints the energy spectrum for a harmonic oscillator
    with a quartic perturbation, based on the provided self-energy diagram.
    """

    # The problem asks for the final equation representing the energy spectrum.
    # The derivation is as follows:
    # 1. The self-energy diagram leads to a correction of the oscillator frequency.
    #    The original Hamiltonian has a term (1/2)*m*omega_0^2 * x^2.
    #    The correction from the diagram adds a term (u*hbar / (8*m*omega_0)) * x^2.
    # 2. The new frequency omega_new is found from:
    #    (1/2)*m*omega_new^2 = (1/2)*m*omega_0^2 + u*hbar / (8*m*omega_0)
    #    omega_new^2 = omega_0^2 + u*hbar / (4*m^2*omega_0)
    # 3. To first order in u, the new frequency is:
    #    omega_new approx omega_0 + u*hbar / (8*m^2*omega_0^2)
    # 4. The energy spectrum is E_n = (n + 1/2)*hbar*omega_new.
    #    E_n = (n + 1/2)*hbar*omega_0 + (n + 1/2)*u*hbar^2 / (8*m^2*omega_0^2)

    # Define the numerical coefficients to be printed.
    coeff_n = 1
    coeff_const = 2
    coeff_corr = 8

    print("The energy spectrum E_n predicted by the self-energy diagram is given by:")
    print(f"E_n = (n + {coeff_n}/{coeff_const}) * hbar * omega_0 + (n + {coeff_n}/{coeff_const}) * u * hbar**2 / ({coeff_corr} * m**2 * omega_0**2)")
    print("\nThis result can be interpreted as a shift in the oscillator's frequency, where the energy levels remain equally spaced but with a new spacing.")

solve_energy_spectrum()
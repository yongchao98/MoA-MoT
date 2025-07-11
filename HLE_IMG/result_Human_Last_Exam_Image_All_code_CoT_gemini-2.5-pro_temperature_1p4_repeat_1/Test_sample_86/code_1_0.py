def calculate_energy_spectrum():
    """
    This function calculates and explains the energy spectrum of a harmonic oscillator
    with a quartic perturbation, as predicted by the tadpole self-energy diagram.
    """
    print("Step 1: The Dyson equation relates the corrected Green's function G to the unperturbed one G_0 and the self-energy Sigma.")
    print("G(omega)^-1 = G_0(omega)^-1 - Sigma(omega)")
    print("The poles of G(omega) give the new energy spectrum. The pole location omega' is found by solving:")
    print("omega'^2 = omega_0^2 + (hbar / m) * Sigma")
    print("-" * 50)
    
    print("Step 2: Calculate the self-energy Sigma from the tadpole diagram.")
    print("Applying the Feynman rules, the self-energy is found to be:")
    # The calculation involves the vertex factor -i(u/hbar), combinatorial factor 1/2,
    # and the loop value G_0(0) = <x^2> = hbar / (2*m*omega_0).
    # -i*Sigma = (-i*u / (2*hbar)) * (hbar / (2*m*omega_0))
    # Sigma = u / (4 * m * omega_0)
    print("Sigma = u / (4 * m * omega_0)")
    print("-" * 50)

    print("Step 3: Calculate the new frequency omega'.")
    print("Substituting Sigma into the equation for the pole:")
    print("omega'^2 = omega_0^2 + (hbar / m) * (u / (4 * m * omega_0))")
    print("omega'^2 = omega_0^2 + (u * hbar) / (4 * m^2 * omega_0)")
    print("\nApproximating for small perturbation u using sqrt(1+x) ~= 1 + x/2:")
    print("omega' = omega_0 * sqrt(1 + (u * hbar) / (4 * m^2 * omega_0^3))")
    print("omega' ~= omega_0 * (1 + (u * hbar) / (8 * m^2 * omega_0^3))")
    print("omega' ~= omega_0 + (u * hbar) / (8 * m^2 * omega_0^2)")
    print("-" * 50)

    print("Step 4: Determine the new energy spectrum E_n'.")
    print("The energy levels are given by E_n = hbar * omega' * (n + 1/2).")
    print("Substituting the expression for the corrected frequency omega':")
    print("E_n' ~= hbar * (omega_0 + (u * hbar) / (8 * m^2 * omega_0^2)) * (n + 1/2)")
    print("\nThe final expression for the energy spectrum, up to an overall constant, is:")
    
    # Final equation presented term by term as requested.
    term1 = "hbar * omega_0 * (n + 1/2)"
    term2 = "(u * hbar^2) / (8 * m^2 * omega_0^2)"
    term3 = "(n + 1/2)"
    print(f"E_n' ~= {term1} + {term2} * {term3}")

calculate_energy_spectrum()
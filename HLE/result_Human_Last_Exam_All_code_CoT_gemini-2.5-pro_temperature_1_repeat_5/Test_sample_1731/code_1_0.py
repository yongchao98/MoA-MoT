def solve_bose_equilibrium():
    """
    Calculates and displays the equilibrium values for mean energy and entropy
    for a Bose gas of light quanta (photons) using principles from
    statistical mechanics.
    """

    print("### Plan ###")
    print("This script outlines the derivation of the equilibrium mean energy (E) and entropy (S) for a photon gas.")
    print("The method is based on the principles of statistical mechanics, which find the most probable macroscopic state by maximizing entropy.")
    print("\nSteps:")
    print("1. Maximize the Bose-Einstein entropy for a system with a fixed total energy E and non-conserved particle number (chemical potential mu=0).")
    print("2. The result is the Bose-Einstein distribution for photons: n(e) = 1 / (exp(e/(k_B*T)) - 1).")
    print("3. Convert sums over discrete energy states to integrals using the density of states g(e) for photons in a 3D volume V.")
    print("4. Solve these integrals to find the final expressions for E and S in terms of temperature (T), volume (V), and fundamental constants (k_B, h, c).")
    print("-" * 60)

    print("\n### Derivation and Results ###\n")

    # --- Equilibrium Mean Energy (E) ---
    print("1. Equilibrium Mean Energy (E)")
    print("The total energy E is found by integrating the energy 'e' weighted by the Bose-Einstein distribution over the density of states g(e).")
    print("   E = Integral[0 to infinity] of (e * g(e) / (exp(e / (k_B * T)) - 1)) de")
    print("where the density of states for photons is g(e) = (8 * pi * V / (h^3 * c^3)) * e^2.")

    print("\nSolving this integral yields the Stefan-Boltzmann law for the total energy of a photon gas.")
    print("\n--- Final Equation for Energy ---")
    # Constants and exponents for the Energy equation
    coeff_E_num = 8
    pow_pi_E = 5
    pow_k_B_E = 4
    pow_T_E = 4
    coeff_E_den = 15
    pow_h_E = 3
    pow_c_E = 3

    print(f"E = ({coeff_E_num} * pi^{pow_pi_E} * k_B^{pow_k_B_E} * V * T^{pow_T_E}) / ({coeff_E_den} * h^{pow_h_E} * c^{pow_c_E})")
    print("\nThis is the equilibrium mean energy of the photon gas, where V is volume, T is temperature,")
    print("k_B is the Boltzmann constant, h is the Planck constant, and c is the speed of light.")
    print("-" * 60)

    # --- Equilibrium Entropy (S) ---
    print("\n2. Equilibrium Entropy (S)")
    print("The entropy (S) for a photon gas in equilibrium is found using the thermodynamic relation S = (4/3) * (E / T).")
    print("We substitute the expression for E found above to get the entropy in terms of the system's parameters.")

    print("\n--- Final Equation for Entropy ---")
    # Constants and exponents for the Entropy equation
    coeff_S_num = 32  # This is 4/3 * 8
    pow_pi_S = 5
    pow_k_B_S = 4
    pow_T_S = 3       # This is T^4 / T
    coeff_S_den = 45  # This is 3 * 15
    pow_h_S = 3
    pow_c_S = 3

    print(f"S = ({coeff_S_num} * pi^{pow_pi_S} * k_B^{pow_k_B_S} * V * T^{pow_T_S}) / ({coeff_S_den} * h^{pow_h_S} * c^{pow_c_S})")
    print("\nThis is the equilibrium entropy of the photon gas.")
    print("-" * 60)

if __name__ == '__main__':
    solve_bose_equilibrium()
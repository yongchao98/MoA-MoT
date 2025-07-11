def print_equilibrium_formulas():
    """
    This function prints the equilibrium formulas for the mean energy and entropy
    of a photon gas (Bose case) based on statistical mechanics.
    """
    print("This script provides the equilibrium values (formulas) for mean energy (E) and entropy (S) for the Bose case of light quanta (a photon gas).")
    print("\nThe derivation is based on maximizing the Bose-Einstein entropy subject to energy conservation.")
    print("\n1. Equilibrium Distribution (Planck Distribution)")
    print("The first step is to find the average number of photons, <n_i>, in each energy state 'i' with energy ε_i.")
    print("The result is the Planck distribution (a special case of the Bose-Einstein distribution with zero chemical potential):")
    print("\n  <n_i> = 1 / (exp(ε_i / (k_B * T)) - 1)\n")
    print("Where:")
    print("  k_B : Boltzmann constant")
    print("  T   : Absolute temperature")
    print("  ε_i : Energy of the i-th state")

    print("\n2. Equilibrium Mean Energy (E)")
    print("The total mean energy E is the sum of the energy of each state multiplied by its average occupation number.")
    print("  E = Σ_i <n_i> * ε_i")
    print("\nSubstituting the Planck distribution, the final equation for mean energy is:")
    # Printing each component of the equation
    print("E", "=", "Σ_i", "[", "ε_i", "/", "(", "exp(", "ε_i", "/", "(", "k_B", "*", "T", ")", ")", "-", "1", ")", "]")

    print("\n3. Equilibrium Entropy (S)")
    print("The entropy S can be derived from the thermodynamic relation S = (E - F) / T, where F is the Helmholtz free energy.")
    print("The free energy for this system is F = k_B * T * Σ_i ln(1 - exp(-ε_i / (k_B * T))).")
    print("\nThis leads to the final equation for entropy:")
    # Printing each component of the equation
    print("S", "=", "E", "/", "T", "-", "k_B", "*", "Σ_i", "ln(", "1", "-", "exp(", "-", "ε_i", "/", "(", "k_B", "*", "T", ")", ")", ")")
    print("\nNote: To find numerical values, the discrete sum (Σ_i) is typically replaced by an integral over a continuous density of states g(ε), which depends on the physical volume of the system.")

print_equilibrium_formulas()
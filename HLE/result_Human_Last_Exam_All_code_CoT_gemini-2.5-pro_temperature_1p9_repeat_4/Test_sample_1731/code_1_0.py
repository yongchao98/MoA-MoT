def solve_photon_gas_equilibrium():
    """
    Derives and prints the equilibrium values for mean energy and entropy
    of a photon gas based on principles of statistical mechanics.
    """

    print("--- Equilibrium State of a Photon Gas ---")
    print("Based on large deviation principles, the equilibrium state of a system corresponds to the state of maximum entropy.")
    print("For a gas of photons (bosons with zero chemical potential), this leads to the Planck distribution.")
    print("\nStep 1: The Planck distribution for the mean occupation number of a state with energy ε at temperature T is:")
    print("<n(ε)> = 1 / (exp(ε / (k_B * T)) - 1)")
    print("where k_B is the Boltzmann constant.\n")

    print("Step 2: To find the mean energy density (U/V), we integrate ε * <n(ε)> over the density of states g(ε).")
    print("The density of states for photons is g(ε) = (8 * π * ε^2) / (h^3 * c^3), where h is Planck's constant and c is the speed of light.")
    print("u = U/V = ∫ [ε * <n(ε)> * g(ε)] dε from 0 to infinity")
    print("Performing this integration yields the Stefan-Boltzmann law for energy density.\n")
    
    print("Step 3: To find the entropy density (S/V), we use the thermodynamic relation for a photon gas: S = (4/3) * U / T.")

    print("\n--- Final Equilibrium Values ---")

    print("\nEquilibrium Mean Energy Density (u = U/V):")
    # Using explicit numbers as requested by the prompt
    print("The final equation for the mean energy density as a function of temperature T is:")
    print("u(T) = (8 * π^5 * k_B^4 / (15 * h^3 * c^3)) * T^4")
    
    print("\nTo show each number in the final equation:")
    print("Numerator constant part: 8 * π^5 * k_B^4")
    print("Denominator constant part: 15 * h^3 * c^3")
    print("The energy is proportional to the temperature to the power of: 4")
    
    print("\nEquilibrium Entropy Density (s = S/V):")
    print("The final equation for the entropy density as a function of temperature T is:")
    print("s(T) = (32 * π^5 * k_B^4 / (45 * h^3 * c^3)) * T^3")

    print("\nTo show each number in the final equation:")
    print("Numerator constant part: 32 * π^5 * k_B^4")
    print("Denominator constant part: 45 * h^3 * c^3")
    print("The entropy is proportional to the temperature to the power of: 3")


if __name__ == '__main__':
    solve_photon_gas_equilibrium()

def print_equilibrium_values():
    """
    This script presents the derived formulas for the equilibrium mean energy (E)
    and entropy (S) for a photon gas (Bose case of light quanta).
    
    The derivation follows from maximizing the system's entropy subject to a
    fixed total energy, a principle underpinned by large deviation theory.
    """
    
    print("The equilibrium values for a photon gas are derived by maximizing entropy.")
    print("This yields the following formulas for the mean energy (E) and entropy (S).")
    print("-" * 60)
    
    # Define symbolic variables for clarity
    constants_desc = (
        "V:   Volume\n"
        "T:   Temperature\n"
        "k_B: Boltzmann constant\n"
        "h:   Planck's constant\n"
        "c:   Speed of light\n"
        "pi:  The constant π"
    )
    print("Constants used in the equations:\n" + constants_desc)
    
    # --- Mean Energy (E) ---
    print("-" * 60)
    print("Equilibrium Mean Energy (E):")
    
    equation_E = "E = (8 * pi^5 * k_B^4 * V * T^4) / (15 * h^3 * c^3)"
    print(f"\n  {equation_E}\n")
    
    print("The integer numbers appearing in the energy equation are:")
    print("  - In the numerator: 8 (from polarization and solid angle), 5 (exponent of pi), 4 (exponent of T and k_B)")
    print("  - In the denominator: 15 (from the integral evaluation), 3 (exponent of h and c)")

    # --- Entropy (S) ---
    print("\n" + "-" * 60)
    print("Equilibrium Entropy (S):")

    equation_S = "S = (32 * pi^5 * k_B^4 * V * T^3) / (45 * h^3 * c^3)"
    print(f"\n  {equation_S}\n")

    print("The integer numbers appearing in the entropy equation are:")
    print("  - In the numerator: 32 (from 4/3 * 8), 5 (exponent of pi), 4 (exponent of k_B), 3 (exponent of T)")
    print("  - In the denominator: 45 (from 3 * 15), 3 (exponent of h and c)")
    print("-" * 60)

# Run the function to display the results
print_equilibrium_values()
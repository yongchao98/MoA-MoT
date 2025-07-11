def absorption_cross_section_equations():
    """
    This script provides the equations for the absorption cross-section of a
    molecular chain interacting with an ultrashort Gaussian laser pulse,
    based on first-order time-dependent perturbation theory.

    It covers two scenarios:
    a) Molecules with no interaction.
    b) Molecules with nearest-neighbor interaction.
    """

    print("=" * 70)
    print("Equation for Absorption Cross-Section of a Molecular Chain")
    print("=" * 70)

    print("\n--- Explanation of Symbols ---\n")
    print("  sigma(omega) : Absorption cross-section as a function of photon frequency.")
    print("  C            : A proportionality constant containing fundamental physical constants.")
    print("  omega        : The angular frequency of the incident light (rad/s).")
    print("  N            : The total number of molecules in the chain.")
    print("  E_0          : The excitation energy of a single, isolated molecule (in Joules).")
    print("  mu_eg        : The magnitude of the transition dipole moment for a single molecule.")
    print("  tau          : The duration of the Gaussian-shaped laser pulse (in seconds).")
    print("  J            : The nearest-neighbor interaction energy, or coupling constant (in Joules).")
    print("  hbar         : The reduced Planck constant.")
    print("  k            : The quantum number for the exciton states (k = 1, 2, ..., N).")
    print("  n            : The position index for a molecule in the chain (n = 1, 2, ..., N).")
    print("\n" + "-" * 70)

    # --- Case a) No Interaction ---
    print("\na) Case with No Interaction (J = 0)\n")
    print("In this model, molecules do not interact. The total absorption is the sum of N identical")
    print("single-molecule absorptions, resulting in a single peak at the monomer energy E_0.")

    print("\nEquation for Absorption Cross-Section (sigma_a):")
    print("sigma_a(omega) = C * omega * (Total_Oscillator_Strength) * (Lineshape_Function)\n")

    print("Where each component is:\n")
    print("  Total Oscillator Strength = N * |mu_eg|^2")
    print("  Lineshape Function        = exp( -(omega - E_0/hbar)^2 * tau^2 )")

    print("\nFinal Equation for Case (a):")
    print("\nsigma_a(omega) = C * omega * N * |mu_eg|^2 * exp(-(omega - E_0/hbar)^2 * tau^2)\n")

    print("-" * 70)

    # --- Case b) Nearest-Neighbor Interaction ---
    print("\nb) Case with Nearest-Neighbor Interaction (J != 0)\n")
    print("Here, interactions delocalize the excitation into exciton states with a band of energies.")
    print("Optical selection rules determine which transitions are allowed ('bright' states).")
    print("The cross-section is a sum over all possible transitions to these exciton states.")

    print("\nEquation for Absorption Cross-Section (sigma_b):")
    print("sigma_b(omega) = C * omega * Sum_{k=1 to N} [ (Oscillator_Strength_k) * (Lineshape_Function_k) ]\n")

    print("Where the components for each exciton state 'k' are:\n")
    print("1. Exciton Energies (E_k):")
    print("   E_k = E_0 + 2 * J * cos( (pi * k) / (N + 1) )")

    print("\n2. Oscillator Strength (S_k):")
    print("   S_k = |mu_eg|^2 * (2 / (N + 1)) * [ Sum_{n=1 to N} sin( (pi * k * n) / (N + 1) ) ]^2")
    print("   (Note: This is non-zero only for odd values of k, forming the selection rule.)")

    print("\n3. Lineshape Function (L_k):")
    print("   L_k(omega) = exp( -(omega - E_k/hbar)^2 * tau^2 )")

    print("\nFinal Equation for Case (b):")
    print("\nsigma_b(omega) = C * omega * Sum_{k=1 to N} { S_k * exp(-(omega - E_k/hbar)^2 * tau^2) }")
    print("   (with E_k and S_k defined as above)")
    print("\n" + "=" * 70)


# Execute the function to print the explanation and equations.
absorption_cross_section_equations()
<<<
def generate_absorption_equations():
    """
    This function formulates and prints the equations for the absorption
    cross-section for a chain of molecules under two conditions.
    """

    # --- Print the explanation of the general formula and terms ---
    print("This program provides the equations for the absorption cross-section \u03C3(E).")
    print("The model is a chain of molecules absorbing an ultrashort Gaussian laser pulse.")
    print("The framework is the first-order time-dependent perturbation theory.")
    print("\nThe general form of the equation is:")
    print("\u03C3(E) = C * E * |d_fi|^2 * G(E)")
    print("where:")
    print("  E      = Photon energy")
    print("  C      = Proportionality constant including fundamental constants")
    print("  |d_fi|^2 = Squared transition dipole moment from initial to final state")
    print("  G(E)   = Gaussian lineshape function due to the ultrashort pulse duration")
    print("\nSymbol definitions:")
    print("  N      = Number of molecules in the chain")
    print("  \u03BC      = Transition dipole moment of a single molecule")
    print("  E_exc  = Excitation energy of a single molecule")
    print("  \u03C4      = Duration of the Gaussian laser pulse")
    print("  \u0127      = Reduced Planck's constant")
    print("  J      = Near-neighbor interaction (coupling) energy\n")

    # --- Case a) No interaction between molecules ---
    print("---")
    print("a) Case with NO interaction between molecules (J=0):")
    print("The final states are N-fold degenerate with energy E_exc.")
    print("The total transition strength is the incoherent sum: N * \u03BC\u00B2.")

    # Construct and print the equation for case a)
    # The \u00B2 is the superscript '2'.
    equation_a = (
        "  \u03C3_a(E) = C * E * (N * \u03BC\u00B2) * "
        "exp[ -((E_exc - E)\u00B2 * \u03C4\u00B2) / \u0127\u00B2 ]"
    )
    print("\nThe equation for the absorption cross-section is:")
    print(equation_a)
    print("---")


    # --- Case b) Near-neighbor interaction ---
    print("\nb) Case with interaction between NEAR-NEIGHBORS (J\u22600):")
    print("The interaction creates a band of delocalized exciton states.")
    print("A selection rule allows transitions only to the k=0 exciton state.")
    print("The energy of this state is shifted to E_exc + 2*J.")
    print("The transition strength is coherently enhanced (superradiance) to N * \u03BC\u00B2.")

    # Construct and print the equation for case b)
    # The \u2260 is the 'not equal to' symbol.
    equation_b = (
        "  \u03C3_b(E) = C * E * (N * \u03BC\u00B2) * "
        "exp[ -((E_exc + 2*J - E)\u00B2 * \u03C4\u00B2) / \u0127\u00B2 ]"
    )
    print("\nThe equation for the absorption cross-section is:")
    print(equation_b)
    print("---")


# Execute the function to print the results
generate_absorption_equations()
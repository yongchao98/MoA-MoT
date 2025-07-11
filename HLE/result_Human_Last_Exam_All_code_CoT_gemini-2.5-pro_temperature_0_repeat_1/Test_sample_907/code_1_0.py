def display_absorption_cross_section_equations():
    """
    This function prints the equations for the absorption cross-section
    for a chain of molecules under two different conditions.
    """

    # --- Introduction ---
    print("This script provides the equations for the absorption cross-section (\u03c3(\u03c9))")
    print("of a molecular chain interacting with an ultrashort Gaussian laser pulse,")
    print("derived using first-order time-dependent perturbation theory.")
    print("-" * 70)

    # --- Case a) No interaction ---
    print("a) Case where interaction between molecules is neglected:")
    print("\n   The total absorption is the sum of absorptions by N independent molecules.")
    print("   The equation is a single Gaussian function centered at the molecular transition frequency \u03c9_fi.")
    print("\n   Equation:")
    print("   \u03c3_a(\u03c9) = N * [ (2 * \u221a\u03c0 * \u03c4 * \u03c9_fi) / (\u210f * c * \u03b5\u2080) ] * |\u20d7d_fi|\u00b2 * exp(-(\u03c9_fi - \u03c9)\u00b2 * \u03c4\u00b2)")
    print("\n   Where:")
    print("   \u03c3_a(\u03c9) : Absorption cross-section as a function of laser frequency \u03c9.")
    print("   N       : Number of molecules in the chain.")
    print("   \u03c4       : Duration of the Gaussian laser pulse.")
    print("   \u03c9_fi     : Transition frequency of a single molecule, (\u0395_f - \u0395_i) / \u210f.")
    print("   \u20d7d_fi     : Transition dipole moment of the molecule, <f|\u20d7d|i>.")
    print("   \u210f       : Reduced Planck constant.")
    print("   c       : Speed of light.")
    print("   \u03b5\u2080       : Vacuum permittivity.")
    print("   exp()   : The exponential function.")
    print("-" * 70)

    # --- Case b) Near-neighbor interaction ---
    print("b) Case where near-neighbor interaction is considered:")
    print("\n   The interaction delocalizes the excitations into a band of Frenkel exciton states.")
    print("   The absorption spectrum is a sum of Gaussians, each corresponding to a transition")
    print("   from the ground state |G> to an optically allowed exciton state |k>.")
    print("\n   Equation:")
    print("   \u03c3_b(\u03c9) = [ (2 * \u221a\u03c0 * \u03c4) / (\u210f * c * \u03b5\u2080) ] * \u2211_{k=1,3,...}^N [ \u03c9_k * |\u20d7d_Gk|\u00b2 * exp(-(\u03c9_k - \u03c9)\u00b2 * \u03c4\u00b2) ]")
    print("\n   Where:")
    print("   \u03c3_b(\u03c9) : Absorption cross-section.")
    print("   \u2211_{k=1,3,...} : Sum over optically allowed exciton states (typically odd integers for a linear chain).")
    print("   \u03c9_k     : Transition frequency to the k-th exciton state, (\u0395_k - \u0395_G) / \u210f.")
    print("   \u0395_k     : Energy of the k-th exciton state. For a simple chain with coupling J,")
    print("             \u0395_k = \u0395_mol + 2*J*cos(\u03c0*k / (N+1)).")
    print("   \u20d7d_Gk     : Transition dipole moment from the ground state |G> to exciton state |k>.")
    print("             It is calculated as <k| \u2211_n \u20d7d_n |G>.")
    print("-" * 70)

if __name__ == '__main__':
    display_absorption_cross_section_equations()
def get_absorption_equations():
    """
    This function prints the equations for the absorption cross-section
    for a chain of molecules under two different conditions as requested.
    """
    print("--- Absorption Cross-Section Equations ---")
    
    # Case a: Non-interacting molecules
    print("\na) The interaction between molecules can be neglected.")
    print("   -------------------------------------------------")
    print("In this case, the total absorption is the incoherent sum of the absorption of N identical molecules.")
    print("\nEquation:")
    print("  sigma(omega) = C * N * omega_0 * mu^2 * exp(-(omega - omega_0)^2 * tau^2)")
    
    print("\nWhere the symbols represent:")
    print("  sigma(omega): Absorption cross-section as a function of laser frequency")
    print("  C          : A proportionality constant containing fundamental physical constants")
    print("  N          : The total number of molecules in the chain")
    print("  omega_0    : The transition frequency of a single, isolated molecule")
    print("  mu         : The magnitude of the transition dipole moment of a single molecule")
    print("  omega      : The central frequency of the Gaussian laser pulse")
    print("  tau        : A parameter related to the duration of the Gaussian laser pulse (shorter pulse -> smaller tau -> wider spectrum)")
    print("  exp()      : The exponential function, e^x. The exponent has a value of 1.")
    print("\nNumbers in the equation:")
    print("  2 : In the exponent, this number represents the quadratic dependence on frequency difference.")
    
    # Case b: Interacting molecules
    print("\n\nb) The interaction between near-neighbors should be considered.")
    print("   ------------------------------------------------------------")
    print("In this case, molecular interactions lead to delocalized exciton states. For a simple chain, a single 'superradiant' state dominates the absorption.")
    print("\nEquation:")
    print("  sigma(omega) = C * N * (omega_0 + 2*J/hbar) * mu^2 * exp(-(omega - (omega_0 + 2*J/hbar))^2 * tau^2)")
    
    print("\nWhere the symbols represent:")
    print("  sigma(omega), C, N, omega_0, mu, omega, tau, exp() are defined as above.")
    print("  J          : The coupling energy (interaction integral) between adjacent molecules. Can be positive or negative.")
    print("  hbar       : The reduced Planck's constant")
    
    print("\nNumbers in the equation:")
    print("  2 : This number appears twice. Once as a factor to the coupling energy J, arising from interactions with two nearest neighbors in the chain. The second is in the exponent, representing the quadratic dependence.")

if __name__ == '__main__':
    get_absorption_equations()

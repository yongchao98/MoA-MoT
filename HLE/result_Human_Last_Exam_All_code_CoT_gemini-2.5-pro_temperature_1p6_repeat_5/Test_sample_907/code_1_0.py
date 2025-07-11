def print_absorption_equations():
    """
    This function prints the equations for the absorption cross-section
    for a chain of molecules in two different cases, as derived from
    first-order time-dependent perturbation theory for a Gaussian pulse.
    """

    # --- Case a) No interaction between molecules ---
    print("---")
    print("a) The interaction between molecules can be neglected.")
    print("---\n")
    print("The equation for the absorption cross-section sigma_a(omega) is:\n")
    
    # Printing each part of the equation symbolically
    equation_a = "sigma_a(omega) = C * N * omega_0 * |mu|^2 * exp(-(omega - omega_0)^2 * tau^2)"
    print(equation_a)
    
    print("\nWhere:")
    print("  sigma_a(omega): Absorption cross-section as a function of laser frequency")
    print("  C:              Proportionality constant including fundamental constants")
    print("  N:              Total number of molecules in the chain")
    print("  omega:          Center frequency of the laser pulse")
    print("  omega_0:        Transition frequency of a single, isolated molecule")
    print("  mu:             Magnitude of the transition dipole moment for a single molecule")
    print("  tau:            Parameter related to the laser pulse duration")
    print("  exp():          The exponential function\n\n")

    # --- Case b) Nearest-neighbor interaction ---
    print("---")
    print("b) The interaction between near-neighbors should be considered.")
    print("---\n")
    print("The equation for the absorption cross-section sigma_b(omega) is:\n")
    
    # Printing each part of the equation symbolically
    equation_b = "sigma_b(omega) = C * N * (omega_0 + 2*J/hbar) * |mu|^2 * exp(-(omega - (omega_0 + 2*J/hbar))^2 * tau^2)"
    print(equation_b)

    print("\nWhere:")
    print("  sigma_b(omega): Absorption cross-section as a function of laser frequency")
    print("  C, N, omega,")
    print("  omega_0, mu, tau: As defined in case (a)")
    print("  J:              Nearest-neighbor coupling energy (interaction strength)")
    print("  hbar:           Reduced Planck's constant")
    print("  exp():          The exponential function")

if __name__ == '__main__':
    print_absorption_equations()
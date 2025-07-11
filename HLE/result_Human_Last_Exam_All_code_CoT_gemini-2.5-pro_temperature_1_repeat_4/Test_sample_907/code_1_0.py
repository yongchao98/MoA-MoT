def absorption_cross_section_equations():
    """
    This function explains and prints the equations for the absorption
    cross-section of a molecular chain interacting with a Gaussian laser pulse,
    based on first-order time-dependent perturbation theory.
    """

    # Introduction to the general theory
    print("The absorption cross-section, sigma(omega), for transitions induced by a Gaussian laser pulse is considered.")
    print("Within the first-order time-dependent perturbation theory, the general formula is:")
    print("\nsigma(omega) = C * sum_over_f [ omega_f * |mu_f|^2 * exp(-(omega - omega_f)^2 * tau^2) ]\n")
    print("where:")
    print("  C         = (sqrt(pi) * tau) / (hbar * c * epsilon_0), a constant of proportionality.")
    print("  omega     = Central frequency of the laser pulse.")
    print("  tau       = Duration of the Gaussian pulse (related to its width).")
    print("  hbar      = Reduced Planck's constant.")
    print("  c         = Speed of light.")
    print("  epsilon_0 = Vacuum permittivity.")
    print("  The sum is over all possible final excited states |f>.")
    print("  omega_f   = Transition frequency from the ground state to the final state |f>.")
    print("  mu_f      = Transition dipole moment for the transition to state |f>.")
    print("-" * 70)

    # Case a) No interaction
    print("\na) Case: The interaction between molecules can be neglected.\n")
    print("In this scenario, the N molecules in the chain are independent. The total absorption is the sum of N identical transitions.")
    print("Each transition has the frequency of a single molecule, omega_0, and a dipole moment of mu_mol.")
    print("\nThe resulting equation for the absorption cross-section is:\n")
    
    equation_a = "sigma_a(omega) = ( (sqrt(pi) * tau) / (hbar * c * epsilon_0) ) * N * omega_0 * |mu_mol|^2 * exp( -(omega - omega_0)^2 * tau^2 )"
    print(equation_a)
    
    print("\nwhere:")
    print("  N         = Total number of molecules in the chain.")
    print("  omega_0   = The transition frequency of a single, isolated molecule.")
    print("  mu_mol    = The transition dipole moment of a single molecule.")
    print("-" * 70)

    # Case b) Near-neighbor interaction
    print("\nb) Case: The interaction between near-neighbors should be considered.\n")
    print("With near-neighbor coupling (J), the excitations delocalize to form Frenkel exciton states.")
    print("For incident light with a wavelength much larger than the chain, a selection rule applies: only the k=0 exciton state is optically active (absorbs light).")
    print("This state has its energy shifted by the coupling and coherently combines the dipole moments of all molecules.\n")
    print(" - Exciton Transition Frequency: omega_k0 = omega_0 + 2*J/hbar")
    print(" - Exciton Transition Dipole Moment (squared): |mu_k0|^2 = N * |mu_mol|^2\n")
    print("The absorption spectrum is dominated by this single, superradiant transition.\n")
    print("The resulting equation for the absorption cross-section is:\n")
    
    equation_b = "sigma_b(omega) = ( (sqrt(pi) * tau) / (hbar * c * epsilon_0) ) * N * |mu_mol|^2 * (omega_0 + 2*J/hbar) * exp( -(omega - (omega_0 + 2*J/hbar))^2 * tau^2 )"
    print(equation_b)

    print("\nwhere:")
    print("  J         = The near-neighbor coupling energy. Its sign determines if the peak shifts to higher (J>0) or lower (J<0) energy.")

# Execute the function to print the results
absorption_cross_section_equations()
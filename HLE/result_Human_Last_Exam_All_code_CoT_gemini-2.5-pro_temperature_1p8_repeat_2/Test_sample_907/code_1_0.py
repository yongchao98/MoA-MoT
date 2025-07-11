import math

def get_absorption_cross_section_equations():
    """
    This function explains and provides the equations for the absorption
    cross-section for a chain of molecules absorbing an ultrashort
    Gaussian laser pulse, derived using first-order time-dependent
    perturbation theory.
    """
    # --- Introduction and General Principles ---
    print("Based on first-order time-dependent perturbation theory, the absorption cross-section (sigma) for a system")
    print("interacting with an ultrashort Gaussian laser pulse depends on the overlap between the pulse's spectrum and")
    print("the system's transition energies. For a pulse with central frequency w_0 and duration tau, the cross-section")
    print("for a transition from an initial state |i> to a final state |f> has a Gaussian shape:\n")
    print("    sigma(w_0) ~ w_fi * |d_fi|^2 * exp( -(w_fi - w_0)^2 * tau^2 / 2 )\n")
    print("where w_fi is the transition frequency, and d_fi is the transition dipole moment.\n")

    # --- Variable Definitions ---
    print("The variables used in the equations are:")
    print("    sigma(w_0): Absorption cross-section as a function of laser frequency w_0")
    print("    N:          Number of molecules in the chain")
    print("    w_0:        Central frequency of the Gaussian laser pulse")
    print("    w_eg:       Electronic transition frequency of a single, isolated molecule")
    print("    tau:        Temporal duration parameter of the Gaussian laser pulse (e.g., E(t) ~ exp(-(t/tau)^2))")
    print("    |mu_eg|^2:  Square of the transition dipole moment for a single molecule")
    print("    J:          Coupling energy between nearest-neighbor molecules")
    print("    hbar:       Reduced Planck's constant")
    print("    c:          Speed of light")
    print("    epsilon_0:  Vacuum permittivity")
    print("\n----------------------------------------------------------------------------------\n")

    # --- Case (a): No interaction ---
    print("a) The interaction between molecules can be neglected.\n")
    print("In this scenario, all N molecules are identical and absorb light independently.")
    print("The total cross-section of the chain is the incoherent sum of the cross-sections")
    print("of each molecule. The absorption spectrum is a Gaussian peak centered at the")
    print("single-molecule transition frequency w_eg.\n")

    equation_a = (
        "sigma_a(w_0) = N * (sqrt(2*pi) / (2 * hbar * c * epsilon_0)) "
        "* |mu_eg|^2 * tau * (w_eg) "
        "* exp( -((w_eg - w_0)^2 * tau^2) / 2 )"
    )

    print("The final equation is:")
    print(equation_a)
    print("\n----------------------------------------------------------------------------------\n")

    # --- Case (b): Nearest-neighbor interaction ---
    print("b) The interaction between near-neighbors should be considered.\n")
    print("In this case, the molecular excitations couple to form delocalized states known as excitons.")
    print("A selection rule dictates that light primarily excites the k=0 exciton state. This state is")
    print("delocalized over the entire chain, leading to a cooperative absorption effect (superradiance).")
    print("The absorption peak is shifted by the coupling energy J.\n")
    print("The transition energy for the k=0 exciton is: w_ex = w_eg + 2*J/hbar")
    print("The transition dipole moment squared for this exciton is enhanced: |d_ex|^2 = N * |mu_eg|^2\n")

    equation_b = (
        "sigma_b(w_0) = (sqrt(2*pi) / (2 * hbar * c * epsilon_0)) "
        "* (N * |mu_eg|^2) * tau * (w_eg + 2*J/hbar) "
        "* exp( -(((w_eg + 2*J/hbar) - w_0)^2 * tau^2) / 2 )"
    )
    
    print("The final equation for the entire chain is:")
    print(equation_b)
    print("\n----------------------------------------------------------------------------------\n")
    
# Execute the function to print the solution
get_absorption_cross_section_equations()

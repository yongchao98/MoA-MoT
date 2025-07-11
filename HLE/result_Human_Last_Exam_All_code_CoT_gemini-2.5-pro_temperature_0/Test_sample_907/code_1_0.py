def display_absorption_equations():
    """
    This function prints the equations for the absorption cross-section
    for a chain of molecules interacting with a Gaussian laser pulse,
    based on first-order time-dependent perturbation theory.
    """

    print("This script presents the equations for the absorption cross-section sigma(omega).")
    print("-" * 70)
    print("Explanation of the symbols used in the equations:")
    print("  omega       : Angular frequency of the incident laser pulse.")
    print("  tau         : Duration (standard deviation) of the Gaussian laser pulse.")
    print("  N           : Number of molecules in the chain.")
    print("  d_eg        : Transition dipole moment for a single molecule's excitation.")
    print("  omega_eg    : Transition angular frequency for a single molecule.")
    print("  J           : Near-neighbor coupling energy.")
    print("  hbar        : Reduced Planck's constant.")
    print("  c           : Speed of light in vacuum.")
    print("  epsilon_0   : Vacuum permittivity.")
    print("  exp(x)      : The exponential function e^x.")
    print("-" * 70)

    # Case a) No interaction between molecules
    print("\na) The interaction between molecules can be neglected.")
    print("In this case, the total absorption is the sum of the absorptions of N independent molecules.")
    print("The absorption spectrum is a Gaussian peak centered at the single-molecule transition frequency omega_eg.")
    print("The width of the peak is determined by the pulse duration tau (a shorter pulse gives a broader peak).")
    
    # Constructing the equation string for case a
    # We output each part of the equation for clarity as requested.
    prefactor = "N * (2 * sqrt(pi) * tau) / (hbar * c * epsilon_0)"
    energy_term_a = "omega_eg"
    dipole_term = "|d_eg|^2"
    gaussian_term_a = "exp(-(omega_eg - omega)^2 * tau^2)"
    
    equation_a = f"sigma_a(omega) = ({prefactor}) * ({energy_term_a}) * ({dipole_term}) * ({gaussian_term_a})"
    
    print("\nThe equation for the absorption cross-section is:")
    print(equation_a)
    
    print("-" * 70)

    # Case b) Near-neighbor interaction
    print("\nb) The interaction between near-neighbors should be considered.")
    print("The interaction leads to the formation of delocalized exciton states.")
    print("Due to selection rules, only one exciton state (the k=0 state) is optically active.")
    print("This state's transition energy is shifted by the coupling energy J.")
    print("The absorption is a single Gaussian peak centered at this new exciton frequency, omega_exc.")

    # Constructing the equation string for case b
    energy_term_b = "omega_exc"
    gaussian_term_b = "exp(-(omega_exc - omega)^2 * tau^2)"
    
    equation_b = f"sigma_b(omega) = ({prefactor}) * ({energy_term_b}) * ({dipole_term}) * ({gaussian_term_b})"
    
    print("\nThe equation for the absorption cross-section is:")
    print(equation_b)
    
    print("\nwhere the exciton transition frequency 'omega_exc' is related to the single-molecule")
    print("frequency 'omega_eg' and the coupling energy 'J' by:")
    omega_exc_relation = "omega_eg + 2*J / hbar"
    print(f"omega_exc = {omega_exc_relation}")


if __name__ == '__main__':
    display_absorption_equations()
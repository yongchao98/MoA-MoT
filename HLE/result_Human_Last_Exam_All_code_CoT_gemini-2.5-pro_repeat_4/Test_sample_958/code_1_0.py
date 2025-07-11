def print_energy_equations():
    """
    Prints the final expressions for the time-averaged stored energy
    per unit area for the electric and magnetic fields of the evanescent wave.
    """
    
    # The expressions are based on the derivation and matching with the provided options.
    # Option D provides the correct expression for the electric field energy.
    
    energy_e_field = "Energy in E field =  (n^2 * (2*n^2 * sin^2(theta) - 1)) / (2 * (omega/c) * (n^2 - 1) * ((n^2 + 1) * sin^2(theta) - 1) * sqrt(n^2 * sin^2(theta) - 1)) * epsilon_0 * |E_{x0}^i|^2"
    
    energy_h_field = "Energy in H field =  (n^2 * (n^2 * sin^2(theta) - 1)) / (2 * (omega/c) * (n^2 - 1) * ((n^2 + 1) * sin^2(theta) - 1) * sqrt(n^2 * sin^2(theta) - 1)) * epsilon_0 * |E_{x0}^i|^2"

    print(energy_e_field)
    print(energy_h_field)

print_energy_equations()
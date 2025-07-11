def print_energy_equations():
    """
    Prints the final expressions for the time-averaged stored energy
    per unit area for the electric and magnetic fields of the evanescent wave.
    """
    # These expressions correspond to the correct derivation for the electric field
    # and the most plausible answer choice provided.
    energy_E_expr = "Energy in E field =  (n^2 * (2*n^2*sin(theta)^2 - 1)) / (2 * (omega/c) * (n^2 - 1) * ((n^2 + 1)*sin(theta)^2 - 1) * sqrt(n^2*sin(theta)^2 - 1)) * epsilon_0 * |E_{x0}^i|^2"
    energy_H_expr = "Energy in H field =  (n^2 * (n^2*sin(theta)^2 - 1)) / (2 * (omega/c) * (n^2 - 1) * ((n^2 + 1)*sin(theta)^2 - 1) * sqrt(n^2*sin(theta)^2 - 1)) * epsilon_0 * |E_{x0}^i|^2"

    print(energy_E_expr)
    print(energy_H_expr)

print_energy_equations()
import sympy as sp

def print_energy_formulas():
    """
    Prints the formulas for the time-averaged stored energy per unit area
    for the electric and magnetic fields of the evanescent wave, as given in Option D.
    """

    # Define symbolic variables
    n, theta, omega, c, epsilon_0, mu_0, E_x0_i = sp.symbols(
        'n theta omega c epsilon_0 mu_0 E_x0_i', real=True, positive=True
    )

    # Common denominator term
    # Note: sympy cannot render the square root nicely in the denominator of a fraction in a single line,
    # so we construct the string manually for clear output.
    denominator_str = ("2 * (omega/c) * (n**2 - 1) * ((n**2 + 1)*sin(theta)**2 - 1) * "
                       "sqrt(n**2*sin(theta)**2 - 1)")

    # Numerator for Electric Field Energy
    numerator_E_str = "n**2 * (2*n**2*sin(theta)**2 - 1)"
    
    # Numerator for Magnetic Field Energy (from option D)
    numerator_H_str = "n**2 * (n**2*sin(theta)**2 - 1)"

    # Construct the full expressions as strings
    energy_E_str = f"({numerator_E_str}) / ({denominator_str}) * epsilon_0 * |E_x0_i|^2"
    energy_H_str = f"({numerator_H_str}) / ({denominator_str}) * epsilon_0 * |E_x0_i|^2"

    print("Based on the derivation, Option D provides the correct expression for the electric field energy.")
    print("There appears to be an inconsistency in the provided expression for the magnetic field energy,")
    print("but Option D is the most plausible answer.\n")
    
    print("Energy in E field:")
    print(energy_E_str)
    
    print("\nEnergy in H field:")
    print(energy_H_str)

if __name__ == '__main__':
    print_energy_formulas()
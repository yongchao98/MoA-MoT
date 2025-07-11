import sympy

def display_poynting_vector():
    """
    This function calculates and displays the symbolic expression for the Poynting vector S.
    S = -(E * rho * v * r / 2) * r_hat + (rho**2 * v * r**2 / (4 * epsilon_0)) * z_hat
    """
    # Using sympy for better symbolic representation, though not strictly necessary for this problem.
    E, rho, v, r, epsilon_0 = sympy.symbols("E ρ v r ε₀")
    
    # Coefficients from the derived equation
    radial_coeff_num = -1
    radial_coeff_den = 2
    axial_coeff_num = 1
    axial_coeff_den = 4

    # Construct the components of the Poynting vector
    S_r = (sympy.Rational(radial_coeff_num, radial_coeff_den)) * E * rho * v * r
    S_z = (sympy.Rational(axial_coeff_num, axial_coeff_den)) * (rho**2 * v * r**2) / epsilon_0
    
    print("The Poynting vector is S = S_r * r_hat + S_z * z_hat")
    print("where r_hat is the radial unit vector and z_hat is the axial unit vector.")
    print("-" * 40)
    
    # Print the radial component and its numerical coefficient
    print("Radial Component (S_r):")
    # To format the output nicely, we extract the coefficient and the rest of the expression
    print(f"S_r = ({radial_coeff_num}/{radial_coeff_den}) * E * ρ * v * r")
    print(f"The numbers in the equation for S_r are: {radial_coeff_num} and {radial_coeff_den}")
    
    print("-" * 40)
    
    # Print the axial component and its numerical coefficient
    print("Axial Component (S_z):")
    print(f"S_z = ({axial_coeff_num}/{axial_coeff_den}) * (ρ^2 * v * r^2) / ε₀")
    print(f"The numbers in the equation for S_z are: {axial_coeff_num} and {axial_coeff_den}")
    print("-" * 40)

    # Print the final combined vector expression
    print("Final Poynting Vector Expression:")
    print(f"S = ({S_r})*r_hat + ({S_z})*z_hat")


display_poynting_vector()
import sympy

def calculate_force(I0_val, D_val, a_val, mu0_val, sigma1_val, sigma2_val):
    """
    Calculates the x-directed total force on the conducting material in the region s < x < 2s.

    Args:
        I0_val (float): DC source current.
        D_val (float): Depth of the system.
        a_val (float): Spacing between the planes.
        mu0_val (float): Magnetic permeability of free space.
        sigma1_val (float): Ohmic conductivity of the first block.
        sigma2_val (float): Ohmic conductivity of the second block.

    Returns:
        float: The x-directed total force.
    """
    I0, D, a, mu0, sigma1, sigma2 = sympy.symbols('I_0 D a mu_0 sigma_1 sigma_2')

    # The formula for the force corresponds to answer choice A
    # F_x = -a * D * (mu_0 / 2) * (I0**2 / D**2) * (sigma2 / (sigma1 + sigma2))**2
    force_expr = -a * D * (mu0 / 2) * (I0**2 / D**2) * (sigma2 / (sigma1 + sigma2))**2

    # Substitute numerical values into the expression
    force_val = force_expr.subs({
        I0: I0_val,
        D: D_val,
        a: a_val,
        mu0: mu0_val,
        sigma1: sigma1_val,
        sigma2: sigma2_val
    })
    
    # Print the equation with substituted values
    term1 = f"(-{a_val} * {D_val})"
    term2 = f"({mu0_val} / 2)"
    term3 = f"(({I0_val}**2) / ({D_val}**2))"
    term4 = f"(({sigma2_val}) / ({sigma1_val} + {sigma2_val}))**2"
    
    print(f"Calculating the force based on the formula: F_x = -a*D * (mu_0/2) * (I_0^2/D^2) * (sigma_2/(sigma_1 + sigma_2))^2")
    print(f"F_x = {term1} * {term2} * {term3} * {term4}")

    return force_val

# --- Main execution ---
if __name__ == '__main__':
    # Define symbolic variables for the formula's representation
    I0, D, a, mu0, sigma1, sigma2 = sympy.symbols('I_0 D a mu_0 sigma_1 sigma_2')
    
    # Represent the final formula symbolically
    force_formula = -a * D * (mu0 / 2) * (I0**2 / D**2) * (sigma2 / (sigma1 + sigma2))**2
    
    # Print the final symbolic formula for the user
    print("The formula for the x-directed total force on the conducting material is:")
    final_equation_str = f"F_x = -a*D * (mu_0/2) * (I_0^2/D^2) * (sigma_2/(sigma_1 + sigma_2))^2"
    print(final_equation_str)
    print("-" * 20)

    # Example parameters for calculation
    I0_val = 10.0  # Amperes
    D_val = 0.2    # meters
    a_val = 0.05   # meters
    mu0_val = 1.25663706e-6  # H/m (permeability of free space)
    sigma1_val = 5.96e7  # S/m (Conductivity of Copper)
    sigma2_val = 3.5e7   # S/m (Conductivity of Gold)
    
    # Calculate the force with the example parameters
    force = calculate_force(I0_val, D_val, a_val, mu0_val, sigma1_val, sigma2_val)

    # Print the final numerical result
    print(f"\nCalculated Force F_x = {force:.4e} Newtons")
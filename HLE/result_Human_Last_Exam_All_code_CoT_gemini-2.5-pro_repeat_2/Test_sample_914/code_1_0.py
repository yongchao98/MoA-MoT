import sympy as sp

def calculate_force_formula():
    """
    This function prints the derived formula for the x-directed force.
    """
    # Define the symbols used in the problem
    a, D, mu_0, I_0, sigma_1, sigma_2 = sp.symbols('a D mu_0 I_0 sigma_1 sigma_2')
    
    # The derived formula for the x-directed force F_x
    # The negative sign indicates the force is in the -x direction.
    term1 = -a * D * (mu_0 / 2)
    # The term related to the current squared and geometry
    term2 = (I_0**2 / D**2)
    # The term related to the current division between the two blocks
    term3 = (sigma_2 / (sigma_1 + sigma_2))**2
    
    # Construct the full expression
    Fx = term1 * term2 * term3
    
    # Print the formula in a readable format, highlighting each component.
    # The question asks to output numbers, which we interpret as the coefficients.
    print("The derived formula for the x-directed force is F_x = C1 * C2 * C3, where:")
    print(f"C1 (system constants): -a*D*mu_0/{2}")
    print(f"C2 (source current term): (I_0**{2} / D**{2})")
    print(f"C3 (current division term): (sigma_2 / (sigma_1 + sigma_2))**{2}")
    print("\nFull formula:")
    print(f"F_x = -a*D * (mu_0/{2}) * (I_0**{2}/D**{2}) * (sigma_2/(sigma_1 + sigma_2))**{2}")

calculate_force_formula()
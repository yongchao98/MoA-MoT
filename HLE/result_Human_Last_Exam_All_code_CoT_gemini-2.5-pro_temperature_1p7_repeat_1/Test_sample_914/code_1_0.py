import sympy

def display_force_formula():
    """
    Defines the system parameters as symbolic variables and prints the formula for the
    x-directed force on the specified conducting block.
    """
    # Define the symbolic variables for the system parameters
    a = sympy.Symbol('a', positive=True)
    D = sympy.Symbol('D', positive=True)
    I_0 = sympy.Symbol('I_0')
    mu_0 = sympy.Symbol('mu_0', positive=True)
    sigma_1 = sympy.Symbol('sigma_1', positive=True)
    sigma_2 = sympy.Symbol('sigma_2', positive=True)

    # Based on the analysis, Choice A represents the force on the second block.
    # F_x = -aD * (mu_0 / 2) * (I_0^2 / D^2) * (sigma_2 / (sigma_1 + sigma_2))^2
    # Simplified: F_x = -(a * mu_0 * I_0**2) / (2 * D) * (sigma_2 / (sigma_1 + sigma_2))**2

    # We will construct the formula from the answer choice A to print its structure.
    term1 = a * D * mu_0 / 2
    term2 = (I_0 / D)**2
    term3 = (sigma_2 / (sigma_1 + sigma_2))**2
    
    Fx = -term1 * term2 * term3

    print("The x-directed total force on the conducting material in the region s < x < 2s is given by the formula:")
    
    # Print the equation in a structured way to show all components
    # The final print statement is constructed to match choice A's structure precisely
    print(f"Fx = - (a*D*{mu_0}/2) * ({I_0}^2/{D}^2) * ({sigma_2}/({sigma_1} + {sigma_2}))^2")
    
    # Using sympy.pretty_print for a more mathematical representation
    print("\nSymbolic representation:")
    sympy.pprint(sympy.Eq(sympy.Symbol('F_x'), Fx), use_unicode=True)

display_force_formula()
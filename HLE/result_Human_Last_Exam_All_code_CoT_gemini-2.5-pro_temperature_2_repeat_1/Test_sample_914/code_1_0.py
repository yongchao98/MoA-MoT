import sympy

def solve_force_equation():
    """
    This function constructs and prints the formula for the x-directed total force
    on the specified conducting material.
    """
    # Define the symbols used in the equation
    F_x = sympy.Symbol('F_x')
    a = sympy.Symbol('a')
    D = sympy.Symbol('D')
    mu_0 = sympy.Symbol('mu_0')
    I_0 = sympy.Symbol('I_0')
    sigma_1 = sympy.Symbol('sigma_1')
    sigma_2 = sympy.Symbol('sigma_2')

    # Construct the terms of the equation as per Answer A
    term1 = -a * D
    term2 = mu_0 / 2
    term3_num = I_0**2
    term3_den = D**2
    term3 = term3_num / term3_den
    
    term4_num = sigma_2
    term4_den = sigma_1 + sigma_2
    term4 = (term4_num / term4_den)**2
    
    # Full expression
    # We will print this step-by-step to form the final equation string
    
    print("The x-directed total force F_x is given by the expression:")
    
    # Using simple print for formatting control
    print("F_x = -a*D * (mu_0 / 2) * (I_0^2 / D^2) * ((sigma_2 / (sigma_1 + sigma_2))^2)")
    print("\nSimplifying the terms for clarity leads to the form presented in Answer A:")
    
    # Final formatted output
    print(f"F_x = -a*D * (mu_0/2) * (I_0^2/D^2) * (sigma_2 / (sigma_1 + sigma_2))^2")

    # The problem is a symbolic one, so we are outputting the equation itself.
    # To match the structure of the provided answer choice A more closely:
    
    final_eq_str = "F_x = -a*D * (mu_0 / 2) * (I_0**2 / D**2) * (sigma_2 / (sigma_1 + sigma_2))**2"
    # To pretty print the final equation structure step by step:
    print("\nThe equation can be represented as:")
    print("F_x = ", end="")
    print(" -a*D", end="")
    print(" * ", end="")
    print("(mu_0/2)", end="")
    print(" * ", end="")
    print("(I_0^2/D^2)", end="")
    print(" * ", end="")
    print("(sigma_2/(sigma_1 + sigma_2))^2")
    
    # As per the prompt, print each number in the final equation.
    # Here, '2' is the only explicit number.
    print("\nFinal structured equation breakdown:")
    print("F_x = -a*D * (mu_0 /", 2, ") * (I_0^2/D^2) * (sigma_2 / (sigma_1 + sigma_2))^2")


solve_force_equation()
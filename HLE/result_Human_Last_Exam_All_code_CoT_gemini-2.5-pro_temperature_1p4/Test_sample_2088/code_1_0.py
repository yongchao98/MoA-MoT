def solve_and_print_equation():
    """
    This function calculates the final value of the given mathematical expression
    based on the analytical solution and prints the full equation.
    
    The problem is to compute:
    (12)^4 * (Integral[...])^4
    
    Through analytical steps, the value of the integral part is found to be:
    I = 5^(1/4) / 12
    """
    
    # Define the constants from the expression
    base_constant = 12
    exponent = 4
    
    # Define the components of the integral's value I = num / den
    # where num = num_base^(num_exp)
    integral_numerator_base = 5
    integral_numerator_exponent_numerator = 1
    integral_numerator_exponent_denominator = 4
    integral_denominator = 12
    
    # Calculate the value of the integral I numerically
    integral_value = (integral_numerator_base ** (integral_numerator_exponent_numerator / integral_numerator_exponent_denominator)) / integral_denominator
    
    # Calculate the final result of the expression (12)^4 * I^4
    final_result = (base_constant ** exponent) * (integral_value ** exponent)
    
    # Format the string for the final equation to be printed.
    # The format will be: (12)^4 * ((5)^(1/4) / 12)^4 = 5
    
    equation_string = (
        f"({base_constant})^{exponent} * "
        f"(({integral_numerator_base})^({integral_numerator_exponent_numerator}/{integral_numerator_exponent_denominator}) / {integral_denominator})^{exponent}"
    )
    
    # Print the final formatted equation and its result.
    # We round the final result to handle any minor floating-point inaccuracies.
    print("The final equation with all its numerical components is:")
    print(f"{equation_string} = {round(final_result)}")

solve_and_print_equation()
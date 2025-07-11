import sympy

def calculate_stress_ratio():
    """
    Calculates the ratio of maximum tangential stress to internal pressure
    for a thick-walled cylinder where the outer radius is twice the inner radius.
    """
    # Define symbols for inner radius (r_i) and outer radius (r_o)
    r_i, r_o = sympy.symbols('r_i r_o')

    print("The formula for the ratio of maximum tangential stress to internal pressure is:")
    # Using unicode for a nicer display of squares
    print("(r_o² + r_i²) / (r_o² - r_i²)\n")

    # Apply the given condition: outer radius is twice the inner radius
    # r_o = 2 * r_i
    condition_eq = sympy.Eq(r_o, 2 * r_i)
    print(f"Given condition: {condition_eq}\n")

    # The expression for the ratio
    ratio_expr = (r_o**2 + r_i**2) / (r_o**2 - r_i**2)

    # Substitute the condition into the expression
    ratio_with_condition = ratio_expr.subs(r_o, 2 * r_i)
    
    print("Substituting the condition into the formula:")
    print("Ratio = ((2*r_i)² + r_i²) / ((2*r_i)² - r_i²)")
    print("Ratio = (4*r_i² + r_i²) / (4*r_i² - r_i²)")
    print("Ratio = (5*r_i²) / (3*r_i²)\n")
    
    # Simplify the expression to get the final numerical ratio
    final_ratio = sympy.simplify(ratio_with_condition)

    # Extract the numerator and denominator to display the final equation
    numerator, denominator = sympy.fraction(final_ratio)

    print("After canceling the r_i² terms, the final equation for the ratio is:")
    print(f"Ratio = {numerator} / {denominator}")
    
    # Also print the decimal value for clarity
    decimal_value = float(final_ratio)
    print(f"\nThe value of the ratio is approximately {decimal_value:.4f}.")

calculate_stress_ratio()
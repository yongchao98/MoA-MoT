import sympy

def solve_stress_ratio():
    """
    Calculates the ratio of maximum tangential stress to internal pressure
    for a thick-walled cylinder where the outer radius is twice the inner radius.
    """
    # Define symbols for the radii and the ratio formula
    a, b = sympy.symbols('a b')
    
    # The formula for the ratio of max tangential stress to internal pressure
    ratio_formula = (b**2 + a**2) / (b**2 - a**2)
    
    # Apply the condition b = 2a
    # We can substitute b with 2*a
    ratio_with_condition = ratio_formula.subs(b, 2*a)
    
    # Simplify the expression
    simplified_ratio = sympy.simplify(ratio_with_condition)
    
    # To fulfill the requirement of showing the numbers in the equation,
    # let's manually substitute values a=1, b=2 into the equation string.
    inner_r_val = 1
    outer_r_val = 2
    
    numerator = outer_r_val**2 + inner_r_val**2
    denominator = outer_r_val**2 - inner_r_val**2
    
    print("The formula for the ratio of maximum tangential stress to internal pressure is:")
    print("Ratio = (b^2 + a^2) / (b^2 - a^2)")
    print("\nGiven the condition that outer radius 'b' is twice the inner radius 'a' (b = 2a).")
    print("We can substitute a=1 and b=2 for simplicity.")
    print("\nThe final equation with substituted values is:")
    print(f"Ratio = ({outer_r_val}^2 + {inner_r_val}^2) / ({outer_r_val}^2 - {inner_r_val}^2)")
    
    print("\nCalculating the numbers in the equation:")
    print(f"Ratio = ({outer_r_val**2} + {inner_r_val**2}) / ({outer_r_val**2} - {inner_r_val**2})")
    
    # Display the final calculation
    print(f"Ratio = {numerator} / {denominator}")

    # The result from sympy is the exact fraction
    print(f"\nThe exact simplified ratio is: {simplified_ratio}")
    # The decimal value
    print(f"The decimal value is approximately: {float(simplified_ratio):.4f}")


solve_stress_ratio()
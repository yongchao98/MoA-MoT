import sympy

def solve_cylinder_stress_ratio():
    """
    Calculates the ratio of maximum tangential stress to internal pressure
    for a thick-walled cylinder where the outer radius is twice the inner radius.
    """
    # Define symbolic variables for inner radius (a) and outer radius (b)
    a, b = sympy.symbols('a b')

    # The formula for the ratio of maximum tangential stress to internal pressure is (b^2 + a^2) / (b^2 - a^2)
    ratio_formula = (b**2 + a**2) / (b**2 - a**2)

    # Apply the given condition: outer radius is twice the inner radius (b = 2a)
    ratio_with_condition = ratio_formula.subs(b, 2*a)

    # Simplify the expression to get the numerical ratio
    final_ratio = sympy.simplify(ratio_with_condition)
    
    # --- Outputting the explanation and calculation ---
    print("The formula for the ratio of maximum tangential stress to internal pressure is:")
    print("Ratio = (b^2 + a^2) / (b^2 - a^2)\n")

    print("Given the condition that the outer radius 'b' is twice the inner radius 'a', we have b = 2a.")
    print("Substituting b = 2a into the formula:\n")
    
    # We can show the substitution by using representative numbers, like a=1 and b=2
    inner_radius_val = 1
    outer_radius_val = 2
    
    b_sq = outer_radius_val**2
    a_sq = inner_radius_val**2
    numerator = b_sq + a_sq
    denominator = b_sq - a_sq
    
    print(f"Ratio = ( (2a)^2 + a^2 ) / ( (2a)^2 - a^2 )")
    print(f"Ratio = ( 4a^2 + a^2 ) / ( 4a^2 - a^2 )")
    print(f"Ratio = ( 5a^2 ) / ( 3a^2 )")
    print(f"After cancelling a^2, the ratio is 5 / 3.\n")

    print("Let's verify with numbers, setting a=1 and b=2:")
    print(f"Ratio = ({outer_radius_val}^2 + {inner_radius_val}^2) / ({outer_radius_val}^2 - {inner_radius_val}^2)")
    print(f"Ratio = ({b_sq} + {a_sq}) / ({b_sq} - {a_sq})")
    print(f"Ratio = {numerator} / {denominator}")

    print(f"\nThe final ratio is {final_ratio}, which is approximately {float(final_ratio):.4f}.")

solve_cylinder_stress_ratio()
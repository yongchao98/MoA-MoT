import sympy

def solve_cylinder_stress():
    """
    Calculates the ratio of maximum tangential stress to internal pressure
    for a thick-walled cylinder where the outer radius is twice the inner radius.
    """
    # Define symbolic variables for the formula
    P_i, a, b = sympy.symbols('P_i a b')

    # The problem states the outer radius (b) is twice the inner radius (a)
    # b = 2a
    relation_b_a = 2 * a

    # Formula for maximum tangential stress (at r=a)
    # sigma_t_max = P_i * (b**2 + a**2) / (b**2 - a**2)
    # We substitute b = 2a into the formula
    numerator_expr = relation_b_a**2 + a**2
    denominator_expr = relation_b_a**2 - a**2

    # Substitute the value of b=2a and simplify
    numerator_val = numerator_expr.subs(b, relation_b_a).simplify()
    denominator_val = denominator_expr.subs(b, relation_b_a).simplify()
    
    # The ratio is sigma_t_max / P_i
    ratio = numerator_val / denominator_val
    
    # Print the steps of the calculation
    print("The formula for the ratio of maximum tangential stress to internal pressure is:")
    print("Ratio = (b^2 + a^2) / (b^2 - a^2)")
    print("\nGiven the condition that the outer radius is twice the inner radius (b = 2a), we substitute this into the formula:")
    # We manually show the substitution for clarity
    print("Ratio = ((2a)^2 + a^2) / ((2a)^2 - a^2)")
    print("Ratio = (4a^2 + a^2) / (4a^2 - a^2)")
    print(f"Ratio = ({numerator_val}) / ({denominator_val})")
    
    # Print the final simplified ratio
    print(f"\nAfter simplifying by cancelling a^2, the final ratio is: {ratio}")
    print(f"As a decimal, this is approximately: {float(ratio):.3f}")

solve_cylinder_stress()
import sympy

def solve_cylinder_stress_ratio():
    """
    Calculates the ratio of maximum tangential stress to internal pressure
    for a thick-walled cylinder where the outer radius is twice the inner radius.
    """
    # Define symbolic variables
    P_i = sympy.Symbol('P_i')  # Internal pressure
    r_i = sympy.Symbol('r_i')  # Inner radius
    r_o = sympy.Symbol('r_o')  # Outer radius

    # Condition: Outer radius is twice the inner radius
    radius_relation = {r_o: 2 * r_i}
    
    # Lamé's formula for maximum tangential stress (at r=r_i, with P_o=0)
    # sigma_t_max = P_i * (r_o**2 + r_i**2) / (r_o**2 - r_i**2)
    
    # We are interested in the ratio sigma_t_max / P_i
    # ratio = (r_o**2 + r_i**2) / (r_o**2 - r_i**2)
    
    numerator = r_o**2 + r_i**2
    denominator = r_o**2 - r_i**2

    # Substitute the radius relation into the numerator and denominator
    numerator_sub = numerator.subs(radius_relation)
    denominator_sub = denominator.subs(radius_relation)

    # Simplify the expressions
    # Numerator becomes: (2*r_i)**2 + r_i**2 = 4*r_i**2 + r_i**2 = 5*r_i**2
    # Denominator becomes: (2*r_i)**2 - r_i**2 = 4*r_i**2 - r_i**2 = 3*r_i**2
    
    numerator_val = 5
    denominator_val = 3
    
    final_ratio = sympy.Rational(numerator_val, denominator_val)
    
    print("The formula for the ratio of maximum tangential stress (σ_t_max) to internal pressure (P_i) is:")
    print("σ_t_max / P_i = (r_o² + r_i²) / (r_o² - r_i²)\n")

    print("Given the condition that the outer radius is twice the inner radius (r_o = 2 * r_i).")
    print("We can substitute this into the formula.\n")

    print("Let's evaluate the numerator and denominator separately:")
    print(f"Numerator = r_o² + r_i²  =>  (2*r_i)² + r_i² = 4*r_i² + 1*r_i² = {numerator_val}*r_i²")
    print(f"Denominator = r_o² - r_i²  =>  (2*r_i)² - r_i² = 4*r_i² - 1*r_i² = {denominator_val}*r_i²\n")

    print("The ratio is the simplified numerator divided by the simplified denominator:")
    print(f"Ratio = ({numerator_val} * r_i²) / ({denominator_val} * r_i²)")
    print(f"After canceling r_i², the final ratio is {numerator_val} / {denominator_val}.\n")

    print(f"The ratio of maximum tangential stress to internal pressure is {final_ratio}.")
    print(f"As a decimal, this is approximately {float(final_ratio):.3f}.")


solve_cylinder_stress_ratio()
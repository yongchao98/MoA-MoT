import sympy

def solve_cylinder_stress_ratio():
    """
    Calculates the ratio of maximum tangential stress to internal pressure
    for a thick-walled cylinder where the outer radius is twice the inner radius.
    """

    # Define the relationship between outer and inner radii
    # Let the inner radius (r_i) be 1 unit for simplicity.
    r_i = 1
    # The outer radius (r_o) is twice the inner radius.
    r_o = 2 * r_i

    # The formula for the ratio of maximum tangential stress to internal pressure is:
    # Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)

    # Calculate the squares of the radii
    r_o_sq = r_o**2
    r_i_sq = r_i**2

    # Calculate the numerator and denominator of the formula
    numerator = r_o_sq + r_i_sq
    denominator = r_o_sq - r_i_sq

    # Print the equation with the values substituted
    print("The formula for the ratio is (r_o^2 + r_i^2) / (r_o^2 - r_i^2)")
    print(f"Given r_o = {r_o} and r_i = {r_i}:")
    print(f"The equation becomes ({r_o}^2 + {r_i}^2) / ({r_o}^2 - {r_i}^2)")
    print(f"= ({r_o_sq} + {r_i_sq}) / ({r_o_sq} - {r_i_sq})")
    print(f"= {numerator} / {denominator}")
    
    # Use sympy to represent the exact fraction
    ratio_fraction = sympy.Rational(numerator, denominator)
    ratio_decimal = float(ratio_fraction)

    print(f"\nThe calculated ratio is {ratio_fraction}, which is approximately {ratio_decimal:.4f}.")
    print("This corresponds to answer choices E (5/3) and H (1.667).")

solve_cylinder_stress_ratio()
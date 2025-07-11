import math

def solve_triangle_area():
    """
    This function derives and prints the formula for the area of triangle T(t).
    
    The key steps of the derivation are:
    1. The rotation of the hexagon does not affect the area of the inscribed triangle.
       The area is calculated in the hexagon's own reference frame.
    2. The triangle T(t) remains equilateral for all t. Its area can be calculated
       using the formula: Area = (sqrt(3)/4) * L^2, where L is the side length.
    3. The side length squared, L(t)^2, is derived as a function of t.
       - Hexagon radius R = 10, so side length s = 10.
       - At t=0, T(0) is an equilateral triangle with side L(0) = 1.5 * s = 15.
       - The vertices of T(t) move a distance v*t = t from the midpoints.
       - The side length squared is L(t)^2 = L(0)^2 + 3*t^2 = 15^2 + 3*t^2.
    4. Substituting L(t)^2 into the area formula gives the final result.
    """
    
    # Constants from the problem
    R = 10
    s = R  # Side length of the hexagon
    v = 1  # Speed of vertices

    # Initial side length of the triangle at t=0
    L0 = 1.5 * s
    
    # The coefficients for the area formula: Area(t) = (sqrt(3)/4) * (L0^2 + 3*t^2)
    # L0^2 = 225
    # The coefficient for t^2 is 3.
    
    l0_squared = L0**2
    t_squared_coeff = 3
    denominator = 4
    sqrt_value = 3
    
    # --- Output the final answer ---
    
    print("Based on the derivation, the area of the triangle T(t) as a function of time t is:")
    
    # Print the formula in a clean, fractional format as requested.
    # This line explicitly outputs each number in the final equation.
    print(f"Area(t) = (sqrt({sqrt_value}) / {denominator}) * ({int(l0_squared)} + {t_squared_coeff}*t^2)")
    
    # For additional clarity, print the formula with decimal coefficients.
    coeff_constant = l0_squared / denominator
    coeff_t_squared = t_squared_coeff / denominator
    
    print("\nAlternatively, expressed with decimal coefficients:")
    print(f"Area(t) = {coeff_constant}*sqrt({sqrt_value}) + {coeff_t_squared}*sqrt({sqrt_value})*t^2")


solve_triangle_area()

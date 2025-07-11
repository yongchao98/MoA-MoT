import math

def find_triangle_area_function():
    """
    This function derives and prints the formula for the area of triangle T(t).
    """
    
    # Step 1: Define the given constants for the problem.
    # Radius of the circle C. This also defines the side length of the inscribed hexagon H.
    R = 10.0
    # Speed of the vertices of T(t) along the sides of H.
    v = 1.0

    # The problem states the hexagon H rotates, but the area of a figure is
    # invariant under rotation. We can ignore the angular velocity ω.

    # Step 2: Calculate the properties of the initial triangle T(0).
    # The side length 's' of the regular hexagon is equal to the radius 'R'.
    s = R
    
    # The vertices of T(0) are the midpoints of three alternating sides of H.
    # Geometrical analysis shows that these points form an equilateral triangle.
    # The side length of this initial triangle, L_0, is (3/2) * s.
    L_0 = (3.0 / 2.0) * s
    
    # The square of the initial side length is:
    L_sq_0 = L_0**2
    
    # Step 3: Determine the side length of the triangle T(t) for any time t.
    # As the vertices move symmetrically from the midpoints, the triangle T(t)
    # remains equilateral.
    # The squared side length L(t)^2 is given by the formula:
    # L(t)^2 = L_0^2 + 3 * (v*t)^2
    
    # We can determine the numerical coefficients for our equation.
    # The constant term comes from the initial squared side length:
    constant_term = int(L_sq_0)
    
    # The coefficient for the t^2 term is 3 * v^2:
    t_squared_coefficient = int(3 * v**2)

    # Step 4: Express the area of the equilateral triangle T(t).
    # The area of an equilateral triangle with side length L is (sqrt(3)/4) * L^2.
    # So, Area(t) = (sqrt(3)/4) * (L_0^2 + 3 * (v*t)^2).
    
    # We will now print the final formula with the calculated numerical values.
    print("The area of the triangle T(t) as a function of time t is given by the equation:")
    print(f"Area(t) = (sqrt(3)/4) * ({constant_term} + {t_squared_coefficient}*t^2)")

find_triangle_area_function()
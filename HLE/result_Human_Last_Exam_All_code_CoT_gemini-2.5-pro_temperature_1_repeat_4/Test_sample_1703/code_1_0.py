import math

def solve_triangle_area():
    """
    This function derives and prints the formula for the area of triangle T(t).
    
    The derivation steps are as follows:
    1.  The hexagon H is inscribed in a circle of radius R=10, so its side length s is also 10.
    2.  At t=0, the triangle T(0) has vertices at the midpoints of alternating sides of H.
        This forms an equilateral triangle with a side length L0 = 1.5 * s = 15.
        The initial squared side length is L0^2 = 15^2 = 225.
    3.  The vertices of T(t) move along the sides of H with speed v=1. The side length L(t)
        of the triangle changes over time. The squared side length is found to be:
        L(t)^2 = L0^2 + 3 * v^2 * t^2
        L(t)^2 = 225 + 3 * (1)^2 * t^2 = 225 + 3 * t^2.
    4.  The rotation of the hexagon does not affect the area of the triangle.
    5.  The area of an equilateral triangle is Area = (sqrt(3)/4) * L^2.
        Substituting L(t)^2, we get the final formula.
    """
    
    # Constants derived from the problem description and geometric analysis
    
    # The constant part of the squared side length, L0^2.
    # L0 = 1.5 * R = 1.5 * 10 = 15. L0^2 = 225.
    constant_term_L_sq = 225
    
    # The coefficient of the t^2 term in the squared side length.
    # This is derived as 3 * v^2, where v = 1.
    t_squared_coefficient_L_sq = 3
    
    # The area formula for an equilateral triangle is (sqrt(3)/4) * L^2.
    # We apply this to L(t)^2 = 225 + 3*t^2.
    
    # The numbers in the final equation
    numerator_sqrt = 3
    denominator_area = 4
    power_of_t = 2
    
    # Print the final equation for the area of T(t)
    print("The area of the triangle T(t) as a function of time t is given by the equation:")
    print(f"Area(t) = (sqrt({numerator_sqrt}) / {denominator_area}) * ({constant_term_L_sq} + {t_squared_coefficient_L_sq} * t^{power_of_t})")

solve_triangle_area()
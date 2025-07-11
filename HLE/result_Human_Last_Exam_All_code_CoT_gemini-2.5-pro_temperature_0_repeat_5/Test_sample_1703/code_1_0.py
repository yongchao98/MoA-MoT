import math

def solve_triangle_area():
    """
    This function calculates and prints the formula for the area of triangle T(t).
    
    The key steps are:
    1. The rotation of the hexagon does not affect the triangle's area, so it can be ignored.
    2. The problem reduces to finding the area of an equilateral triangle whose side length changes with time.
    3. The initial side length squared (at t=0) is a_0^2 = 225.
    4. The side length squared as a function of time t is a(t)^2 = 225 + 3*t^2.
    5. The area of an equilateral triangle is Area = (sqrt(3)/4) * a^2.
    """

    # The constant part of the side length squared equation, a_0^2
    side_squared_const = 225

    # The coefficient of the t^2 term in the side length squared equation
    t_squared_coeff = 3

    # The denominator in the area formula for an equilateral triangle
    area_denominator = 4

    # Print the final equation for the area as a function of t.
    # The format is chosen to clearly show all the numbers in the final equation.
    print("The area of the triangle T(t) as a function of time t is given by the equation:")
    print(f"Area(t) = (sqrt(3) / {area_denominator}) * ({side_squared_const} + {t_squared_coeff} * t^2)")

solve_triangle_area()
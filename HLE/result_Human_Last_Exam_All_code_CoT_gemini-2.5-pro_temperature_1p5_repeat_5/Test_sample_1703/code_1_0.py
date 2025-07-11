def solve_triangle_area():
    """
    This function prints the formula for the area of triangle T(t).
    The formula is derived from the geometric and kinematic properties described in the problem.
    A(t) = (sqrt(3)/4) * (225 + 3*t^2)
    """

    # The numbers in the final equation as derived from the analysis:
    # A(t) = (sqrt(c1)/c2) * (c3 + c4 * t^2)
    c1 = 3   # The number inside the square root.
    c2 = 4   # The denominator.
    c3 = 225 # The constant term, from the initial side length of the triangle.
    c4 = 3   # The coefficient of the t^2 term.

    # Print the final equation for the area A(t) as a function of time t.
    print("The area of the triangle T(t) as a function of time t is given by the formula:")
    print(f"A(t) = (sqrt({c1}) / {c2}) * ({c3} + {c4}*t^2)")

solve_triangle_area()
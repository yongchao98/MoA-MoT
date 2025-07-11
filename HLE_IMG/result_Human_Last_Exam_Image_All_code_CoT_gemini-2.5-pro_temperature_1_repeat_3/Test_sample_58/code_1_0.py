import sympy

def solve_cylinder_height():
    """
    This function derives and prints the formula for the height of the cylinder.
    """
    # Define the symbolic variables for height, radius, and theta.
    h, r, theta = sympy.symbols('h r theta')
    pi = sympy.pi

    # The derivation is based on unrolling the cylinder's surface.
    # A right-angled triangle is formed on this unrolled surface.
    # The sides of the triangle are:
    # - The height of the cylinder, h.
    # - The arc length between two diametrically opposite points, pi * r.
    # The angle theta from the problem is interpreted as the angle in this triangle
    # such that tan(theta) = h / (pi * r).

    # The resulting equation for the height h is:
    # h = pi * r * tan(theta)

    print("The height 'h' of the cylinder is given by the formula:")
    # We print the formula string by string to satisfy the output format.
    print("h", "=", "pi", "*", "r", "*", "tan(theta)")

solve_cylinder_height()
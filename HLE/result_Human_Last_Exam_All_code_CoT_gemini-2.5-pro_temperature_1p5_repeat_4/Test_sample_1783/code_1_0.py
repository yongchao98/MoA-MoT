import sympy

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2(x^2+y^2) over a cone.
    The cone has a base radius of 3 at z=0 and a height of 2.
    """
    # Define the symbolic variables
    r, theta, z = sympy.symbols('r theta z')
    pi = sympy.pi

    # The function in cylindrical coordinates is f = z^2 * r^2.
    # The integrand is f * r (the Jacobian for cylindrical coordinates).
    integrand = z**2 * r**3

    # Define the limits of integration for the cone.
    # Height of the cone H=2, Radius of the base R=3.
    # The radius r at a given height z is r(z) = R * (1 - z/H) = 3 * (1 - z/2).
    r_limit_upper = 3 - (sympy.S(3)/2) * z
    
    # Set up the integration bounds
    # Order of integration: dr, dtheta, dz
    integration_bounds = [
        (r, 0, r_limit_upper),
        (theta, 0, 2 * pi),
        (z, 0, 2)
    ]

    # Calculate the definite triple integral
    # The integrate function takes the integrand and the bounds in order
    result = sympy.integrate(integrand, *integration_bounds)
    
    # The result is (108*pi)/35. We need to print each number.
    # We can extract the numerator and denominator from the fraction.
    if result.has(pi):
        coefficient = result / pi
        p, q = sympy.fraction(coefficient)
        print(f"The integral of f(x,y,z) = z^2 * (x^2 + y^2) over the specified cone is:")
        print(f"Final equation: ({p} * pi) / {q}")
    else:
        # Fallback for non-pi results
        print(f"The integral evaluates to: {result}")

solve_cone_integral()
<<<108*pi/35>>>
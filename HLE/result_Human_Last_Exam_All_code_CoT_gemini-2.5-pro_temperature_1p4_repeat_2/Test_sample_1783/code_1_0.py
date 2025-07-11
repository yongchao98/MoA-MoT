import sympy

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2*(x^2+y^2) over a cone
    with base radius 3 at z=0 and height 2.
    """
    # Define the symbols for cylindrical coordinates
    r, theta, z = sympy.symbols('r theta z')
    
    # Use sympy's symbolic pi
    pi = sympy.pi

    # In cylindrical coordinates, f(r, theta, z) = z^2 * r^2.
    # The volume element dV is r*dr*d(theta)*dz.
    # The full integrand is f * r = z^2 * r^3.
    integrand = z**2 * r**3

    # Define the integration limits for the cone.
    # z (height) goes from 0 to 2.
    # theta (angle) goes from 0 to 2*pi for a full circle.
    # r (radius) depends on z. It goes from 0 to the edge of the cone,
    # which is described by the line r = 3 * (1 - z/2).
    z_limits = (z, 0, 2)
    theta_limits = (theta, 0, 2 * pi)
    r_limits = (r, 0, 3 * (1 - z/2))

    # Calculate the triple integral.
    # We can integrate in any order since the limits are well-defined.
    integral_result = sympy.integrate(integrand, r_limits, theta_limits, z_limits)

    # The problem asks to output each number in the final equation.
    # The result is of the form (A * pi) / B.
    # We will extract A and B.
    if pi in integral_result.args:
        coefficient = integral_result / pi
        numerator = sympy.numer(coefficient)
        denominator = sympy.denom(coefficient)

        print("The integral is calculated in the form (A * pi) / B.")
        print("Final Equation: Integral = ({} * pi) / {}".format(numerator, denominator))
        print("Value of A:", numerator)
        print("Value of B:", denominator)
    else:
        # Fallback for a result that does not contain pi
        print("The result of the integration is:", integral_result)

solve_cone_integral()
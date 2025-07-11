import sympy

def solve_integral():
    """
    This function calculates the integral of f(x,y,z) = z^2(x^2+y^2) over a cone
    with base radius 3 and height 2, centered at the origin.
    """
    # Define the symbols for our cylindrical coordinates
    z, r, theta = sympy.symbols('z r theta')
    pi = sympy.pi

    # The function in cylindrical coordinates is f = z^2 * r^2
    # The volume element dV is r*dr*dtheta*dz
    # So the full integrand is (z^2 * r^2) * r
    integrand = z**2 * r**3

    # Define the limits of integration for the cone.
    # Height H=2, Radius R=3.
    # The radius of the cone's cross-section at height z is r_z = R * (H-z) / H
    r_limit_upper = sympy.S(3)/2 * (2 - z)

    # We set up the integration order as dr, then dtheta, then dz.
    # Integrate with respect to r
    integral_r = sympy.integrate(integrand, (r, 0, r_limit_upper))

    # Integrate with respect to theta
    integral_theta = sympy.integrate(integral_r, (theta, 0, 2 * pi))

    # Integrate with respect to z to get the final result
    final_result = sympy.integrate(integral_theta, (z, 0, 2))

    # To output each number in the final equation, we extract the numerator and
    # denominator of the fraction part of the answer.
    # The result is of the form C * pi. We want to show it as (A * pi) / B
    if final_result.is_rational:
        # The result is a number without pi
        p = final_result
        q = 1
        print(f"The final result is an equation: {p} / {q} = {final_result.evalf():.4f}")
    elif isinstance(final_result, sympy.Mul) and pi in final_result.args:
        # The result is of the form C * pi
        coeff = final_result / pi
        p, q = sympy.fraction(coeff)
        print(f"The final equation for the integral is: ({p} * pi) / {q} = {final_result.evalf():.4f}")
    else:
        # For any other form, just print it
        print(f"The final result of the integral is: {final_result} = {final_result.evalf():.4f}")


solve_integral()
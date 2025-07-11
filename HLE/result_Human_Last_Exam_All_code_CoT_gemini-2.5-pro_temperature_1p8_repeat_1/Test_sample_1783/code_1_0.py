import sympy

def solve_integral():
    """
    This function calculates the integral of f(x,y,z) = z^2*(x^2+y^2) over a cone.
    The cone has a base radius of 3 at z=0 and an apex at (0,0,2).
    """
    # Define the symbolic variables
    r, theta, z = sympy.symbols('r theta z')

    # Define the function in cylindrical coordinates: f = z^2 * r^2
    # The integrand includes the Jacobian 'r', so it becomes z^2 * r^3
    integrand = z**2 * r**3

    # Define the limits of integration for the cone
    z_min, z_max = 0, 2
    r_min = 0
    # The radius 'r' depends on 'z'. r = (3/2)*(2-z)
    r_max = sympy.Rational(3, 2) * (2 - z)
    theta_min, theta_max = 0, 2 * sympy.pi

    # Integrate with respect to r first
    integral_r = sympy.integrate(integrand, (r, r_min, r_max))
    
    # Next, integrate the result with respect to z
    integral_z = sympy.integrate(integral_r, (z, z_min, z_max))
    
    # Finally, integrate the result with respect to theta
    final_integral_value = sympy.integrate(integral_z, (theta, theta_min, theta_max))

    # The result is of the form (numerator * pi) / denominator.
    # We extract these parts to clearly display the final equation.
    if final_integral_value.has(sympy.pi):
        # Extract the coefficient of pi
        coeff = final_integral_value / sympy.pi
        numerator, denominator = sympy.fraction(coeff)
        print(f"The integral is ({numerator}*pi)/{denominator}")
    else:
        # Handle cases where pi might cancel out or not be present
        numerator, denominator = sympy.fraction(final_integral_value)
        if denominator == 1:
            print(f"The integral is {numerator}")
        else:
            print(f"The integral is {numerator}/{denominator}")

solve_integral()
import sympy

def solve_cone_integral():
    """
    This function calculates the integral of f(x,y,z) = z^2*(x^2+y^2)
    over a cone with base radius 3 and height 2, centered at the origin.
    """
    # 1. Define the symbolic variables
    r, theta, z = sympy.symbols('r theta z')
    pi = sympy.pi

    # 2. Define the function and the integrand in cylindrical coordinates
    # The function is f = z^2 * r^2
    # The volume element gives an extra factor of r, so the integrand is z^2 * r^3
    integrand = z**2 * r**3

    # 3. Define the limits of integration for the cone
    # z goes from 0 to height H=2
    # theta goes from 0 to 2*pi for a full cone
    # r goes from 0 to the edge of the cone, which is r = 3*(1 - z/2)
    r_limits = (r, 0, 3 * (1 - z/2))
    theta_limits = (theta, 0, 2 * pi)
    z_limits = (z, 0, 2)

    # 4. Compute the triple integral step-by-step (inside-out)
    integral_over_r = sympy.integrate(integrand, r_limits)
    integral_over_theta = sympy.integrate(integral_over_r, theta_limits)
    final_integral = sympy.integrate(integral_over_theta, z_limits)

    # 5. Format and print the final equation as requested
    # The result is of the form (numerator * pi) / denominator
    numerator, denominator = sympy.fraction(final_integral / pi)
    
    print("The final equation for the integral is:")
    print(f"Integral = ({numerator} * \u03C0) / {denominator}")

solve_cone_integral()
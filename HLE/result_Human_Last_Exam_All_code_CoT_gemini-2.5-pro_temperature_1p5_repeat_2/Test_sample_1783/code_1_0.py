import sympy
from sympy import integrate, pi, symbols

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2(x^2+y^2) over a cone.
    The cone has a base radius of 3 at z=0 and a height of 2.
    """
    # Define the symbols for our integration
    r, theta, z = symbols('r theta z')

    # The function f in cylindrical coordinates is z^2 * r^2.
    # The volume element dV is r*dr*dtheta*dz.
    # The integrand is f * r.
    integrand = z**2 * r**3

    # Define the limits of integration
    # r: 0 to 3 - (3/2)*z
    # theta: 0 to 2*pi
    # z: 0 to 2
    r_limit_upper = 3 - (3 * z) / 2
    theta_limit_upper = 2 * pi
    z_limit_upper = 2
    
    # We will print the setup of the final calculation
    # First, integrate with respect to r
    integral_r = integrate(integrand, (r, 0, r_limit_upper))
    
    # Next, integrate the result with respect to theta
    integral_theta = integrate(integral_r, (theta, 0, theta_limit_upper))

    # Finally, integrate the result with respect to z to get the final answer
    final_integral_result = integrate(integral_theta, (z, 0, z_limit_upper))

    # Deconstruct the final answer to show its components as requested
    # The result is of the form C * pi
    if final_integral_result.is_rational:
        num = final_integral_result.p
        den = final_integral_result.q
        print(f"The final equation is: ({num} * pi) / {den}")
    elif (final_integral_result / pi).is_rational:
        fraction = final_integral_result / pi
        num = fraction.p
        den = fraction.q
        print(f"The final equation is: ({num} * pi) / {den}")
    else:
        # Fallback for unexpected formats
        print(f"The final result is: {final_integral_result}")

    print(f"The numerical value is approximately: {final_integral_result.evalf()}")


solve_cone_integral()
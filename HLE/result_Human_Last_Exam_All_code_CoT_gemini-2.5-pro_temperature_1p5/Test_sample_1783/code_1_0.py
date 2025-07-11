import numpy as np
from scipy import integrate

def solve_integral():
    """
    Calculates the integral of f(x,y,z) = z^2(x^2+y^2) over a cone.
    The cone has a base radius of 3 at z=0 and a height of 2 (vertex at z=2).
    The calculation is performed using triple integration in cylindrical coordinates.
    """

    # The function to integrate in cylindrical coordinates is z^2 * r^2.
    # The volume element dV is r*dr*d(theta)*dz.
    # So, the final integrand is (z^2 * r^2) * r = z^2 * r^3.
    # For scipy's tplquad, the function arguments are in the order of integration:
    # outer, middle, inner. Our order is z, theta, r.
    # So we define the function as f(z, theta, r).
    integrand = lambda r, theta, z: z**2 * r**3

    # The limits for z are from 0 to 2.
    z_min, z_max = 0, 2
    
    # The limits for theta are from 0 to 2*pi.
    # These are functions of the outer variable, z.
    theta_min = lambda z: 0
    theta_max = lambda z: 2 * np.pi

    # The limits for r depend on z: r = 3 * (1 - z/2).
    # These are functions of the outer two variables, z and theta.
    r_min = lambda z, theta: 0
    r_max = lambda z, theta: 3 * (1 - z / 2.0)
    
    print("Calculating the integral of f(x,y,z) = z^2 * (x^2 + y^2) over a cone.")
    print("Using cylindrical coordinates, the integral is:")
    print("I = integral from z=0 to 2 of ( integral from theta=0 to 2*pi of ( integral from r=0 to 3*(1-z/2) of (z^2 * r^3) dr ) d(theta) ) dz")

    # Perform the triple integration.
    # Note: tplquad requires the function signature to match the integration order (z, theta, r).
    # We use a lambda function to reorder the arguments passed to our integrand.
    result, error = integrate.tplquad(
        lambda z, theta, r: integrand(r, theta, z),
        z_min, z_max,
        theta_min, theta_max,
        r_min, r_max
    )

    analytical_result_val = 108 * np.pi / 35
    
    print("\nThe analytical result is the fraction 108*pi/35.")
    print("The final equation and its numerical value are:")
    # As requested, outputting each number in the final equation.
    print("108 * {} / 35 = {}".format(np.pi, result))
    
    # Returning the numerical value for the final answer tag.
    return result

if __name__ == '__main__':
    final_value = solve_integral()
    # print(f"\nFinal numerical answer: {final_value}")
    # print(f'<<<{final_value:.3f}>>>')

solve_integral()
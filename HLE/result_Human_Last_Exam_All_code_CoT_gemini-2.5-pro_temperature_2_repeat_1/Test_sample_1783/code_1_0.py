import sympy

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2(x^2+y^2) over a cone.
    The cone has a base radius of 3 at z=0 and a height of 2.
    """
    # Define the symbols for cylindrical coordinates
    r, theta, z = sympy.symbols('r theta z')

    # Define the function in cylindrical coordinates
    # f(x,y,z) = z^2(x^2+y^2) becomes f(r,theta,z) = z^2 * r^2
    f_cylindrical = z**2 * r**2

    # The volume element in cylindrical coordinates is r*dr*dtheta*dz
    integrand = f_cylindrical * r

    # Define the limits of integration
    # z (height): from 0 to 2
    # theta (angle): from 0 to 2*pi
    # r (radius): from 0 to 3 - (3/2)*z
    r_limit_upper = 3 - (3 * z / 2)
    
    # Calculate the definite triple integral
    # The order of integration is dr, then dtheta, then dz
    integral_value = sympy.integrate(
        integrand,
        (r, 0, r_limit_upper),
        (theta, 0, 2 * sympy.pi),
        (z, 0, 2)
    )

    # To fulfill the request "output each number in the final equation!",
    # we will extract the numerator and denominator of the result.
    num, den = integral_value.as_numer_denom()
    
    # The numerator is of the form C*pi. We extract C.
    coeff = num.as_coeff_Mul()[0]

    print("The final integral is:")
    # Print the equation representing the result
    print(f"({coeff} * \u03C0) / {den}")
    print(f"\nThis evaluates to approximately: {integral_value.evalf()}")

solve_cone_integral()
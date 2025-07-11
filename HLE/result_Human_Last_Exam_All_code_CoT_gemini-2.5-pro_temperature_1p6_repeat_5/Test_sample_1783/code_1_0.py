import sympy as sp

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2(x^2+y^2) over a cone.
    """
    # Define the symbols for cylindrical coordinates
    r, theta, z = sp.symbols('r theta z')

    # The function f(x,y,z) = z^2(x^2+y^2) in cylindrical coordinates is z^2*r^2
    # The integrand is f * r (from the volume element r*dr*dz*d_theta)
    integrand = z**2 * r**3

    # Define the limits of integration for the cone
    # Height H=2, Radius R=3
    # theta: 0 to 2*pi
    # z: 0 to H=2
    # r: 0 to R*(1 - z/H) = 3*(1-z/2) = 3 - 3*z/2
    limits = [
        (r, 0, 3 - (3*sp.S(1)/2)*z),
        (z, 0, 2),
        (theta, 0, 2*sp.pi)
    ]

    # Create the integral expression for display purposes
    integral_expression = sp.Integral(integrand, *limits)

    print("The integral to be solved is:")
    sp.pprint(integral_expression, use_unicode=True)
    print("\n")

    # Evaluate the integral
    result = sp.integrate(integrand, *limits)
    
    # Extract numbers for the final print statement
    num, den = result.as_numer_denom() # e.g., 108*pi, 35
    coeff = num.as_coeff_Mul()[0] # e.g., 108
    
    print(f"The final result of the integration is: {result}")
    
    # "output each number in the final equation!"
    # The final equation is: Integral_Value = (108 * pi) / 35
    # The numbers are 108 and 35.
    print(f"The final equation consists of the numbers {int(coeff)} and {int(den)}, multiplied by pi.")
    
    print(f"The numerical value is approximately: {result.evalf()}")

solve_cone_integral()
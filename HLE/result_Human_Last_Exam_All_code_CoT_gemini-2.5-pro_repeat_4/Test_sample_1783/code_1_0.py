import sympy

def solve_integral():
    """
    Calculates the integral of f(x,y,z) = z^2 * (x^2 + y^2) over a cone
    with base radius 3 and height 2.
    """
    # Define the symbolic variables
    r, theta, z = sympy.symbols('r theta z')
    pi = sympy.pi

    # Define the function in cylindrical coordinates and the volume element
    # f(r, theta, z) = z^2 * r^2
    # Integrand = f * r = z^2 * r^3
    integrand = z**2 * r**3

    # Define the integration limits
    # Height H=2, Radius R=3
    # Upper limit for r: r(z) = R * (H - z) / H = 3 * (2 - z) / 2
    r_upper_limit = sympy.S(3)/2 * (2 - z)
    
    # Perform the triple integration
    # 1. Integrate with respect to r
    integral_r = sympy.integrate(integrand, (r, 0, r_upper_limit))
    
    # 2. Integrate the result with respect to theta
    integral_theta = sympy.integrate(integral_r, (theta, 0, 2 * pi))
    
    # 3. Integrate the result with respect to z
    final_result = sympy.integrate(integral_theta, (z, 0, 2))

    # Extract the components of the final answer for printing
    if final_result.is_rational:
        num = final_result
        den = 1
        has_pi = False
    else:
        num, den = sympy.fraction(final_result / pi)
        has_pi = True

    # Print the final equation as requested
    print("The result of the integral is given by the following equation:")
    if has_pi:
        print(f"({num} * pi) / {den}")
        print("\nBreaking down the final equation:")
        print(f"The numerator is: {num}")
        print(f"The denominator is: {den}")
        print("The equation includes the symbol: pi")
    else:
        print(f"{num} / {den}")
        print("\nBreaking down the final equation:")
        print(f"The numerator is: {num}")
        print(f"The denominator is: {den}")

    print(f"\nThe numerical value is approximately: {final_result.evalf()}")


solve_integral()
<<<108*pi/35>>>
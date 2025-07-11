import sympy

def calculate_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2*(x^2+y^2) over a specified cone
    using symbolic integration in cylindrical coordinates.
    """
    # Define the symbolic variables
    r, theta, z = sympy.symbols('r theta z')

    # Define the function and the volume element product in cylindrical coordinates
    # f(r, theta, z) = z^2 * r^2
    # dV = r * dr * dtheta * dz
    # Integrand = f * r = z^2 * r^3
    integrand = z**2 * r**3

    # Define the limits of integration
    # Height of the cone: H = 2
    # Radius of the base: R = 3
    # The radius of the cone at height z is r(z) = R * (1 - z/H) = 3 * (1 - z/2)
    r_limit = sympy.S(3)/2 * (2 - z)
    
    # Perform the integration sequentially
    # 1. Integrate with respect to r
    integral_r = sympy.integrate(integrand, (r, 0, r_limit))
    
    # 2. Integrate the result with respect to theta
    integral_theta = sympy.integrate(integral_r, (theta, 0, 2 * sympy.pi))
    
    # 3. Integrate the result with respect to z
    final_result = sympy.integrate(integral_theta, (z, 0, 2))

    # The final symbolic result is 108*pi/35.
    # We will display the numbers in the final equation as requested.
    if isinstance(final_result, sympy.Mul) and sympy.pi in final_result.args:
        # Extract the coefficient of pi
        coeff = final_result / sympy.pi
        numerator, denominator = sympy.fraction(coeff)
        print(f"The final integral is equal to the equation: ({numerator} * pi) / {denominator}")
    else:
        print(f"The final integral is: {final_result}")

    # Print the numerical value for the final answer format
    numerical_value = final_result.evalf()
    print(f"The numerical value of the integral is: {numerical_value}")
    
calculate_cone_integral()
<<<9.70424696347374>>>
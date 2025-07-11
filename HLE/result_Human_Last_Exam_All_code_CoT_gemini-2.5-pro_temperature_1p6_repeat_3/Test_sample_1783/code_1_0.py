import sympy

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2*(x^2+y^2) over a cone
    with base radius 3 and height 2, using symbolic integration.
    """
    # Define the symbols for cylindrical coordinates
    r, z, theta = sympy.symbols('r z theta')

    # Define the cone's dimensions
    H = 2  # Height
    R = 3  # Radius

    # The function in cylindrical coordinates is z^2 * r^2
    # The volume element dV is r*dr*dz*dtheta
    # The full integrand is (z^2 * r^2) * r
    integrand = z**2 * r**3

    # Define the upper limit for the radius 'r' as a function of 'z'
    # r = R * (1 - z/H)
    r_upper_limit = R * (1 - z / H)

    # Perform the triple integration step-by-step
    # 1. Integrate with respect to r
    integral_r = sympy.integrate(integrand, (r, 0, r_upper_limit))
    
    # 2. Integrate the result with respect to z
    integral_z = sympy.integrate(integral_r, (z, 0, H))
    
    # 3. Integrate the result with respect to theta
    final_result = sympy.integrate(integral_z, (theta, 0, 2 * sympy.pi))
    
    # Extract the components of the final answer to display them clearly
    # The result is of the form (A/B) * pi
    # We divide by pi to get the fraction part
    fractional_part = final_result / sympy.pi
    num, den = sympy.fraction(fractional_part)

    print(f"The integral of f(x,y,z) = z^2(x^2+y^2) over the given cone is calculated as:")
    print(f"Result = ({num} * pi) / {den}")
    # Also print the numerical value for reference
    print(f"Numerical value: {final_result.evalf()}")

solve_cone_integral()
<<<108*pi/35>>>
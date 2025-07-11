import sympy

def solve_cone_integral():
    """
    This function calculates the integral of f(x,y,z) = z^2(x^2+y^2)
    over the volume of a cone with base radius 3 and height 2.
    """
    # 1. Define the symbols and the function in cylindrical coordinates
    r, theta, z = sympy.symbols('r theta z')
    pi = sympy.pi

    # The function is f = z^2 * r^2. The integrand is f * r (Jacobian).
    integrand = z**2 * r**3

    # 2. Define the integration limits
    # Limit for r: 0 to 3/2 * (2-z)
    r_limit_upper = sympy.Rational(3, 2) * (2 - z)
    # Limit for z: 0 to 2
    z_limit_upper = 2
    # Limit for theta: 0 to 2*pi
    theta_limit_upper = 2 * pi

    # 3. Perform the integration step-by-step and print the results
    print("The integral to be solved is:")
    print(f"∫ from 0 to {theta_limit_upper} d(theta) ∫ from 0 to {z_limit_upper} dz ∫ from 0 to {r_limit_upper} ({integrand}) dr")
    print("-" * 30)

    # Step 1: Integrate with respect to r
    integral_r = sympy.integrate(integrand, (r, 0, r_limit_upper))
    print(f"Step 1: The result of the innermost integral (with respect to r) is:")
    print(f"∫ ({integrand}) dr = {integral_r}")
    print("-" * 30)
    
    # Step 2: Integrate the result with respect to z
    integral_z = sympy.integrate(integral_r, (z, 0, z_limit_upper))
    print(f"Step 2: The result of the middle integral (with respect to z) is:")
    print(f"∫ ({integral_r}) dz = {integral_z}")
    print("-" * 30)

    # Step 3: Integrate the result with respect to theta
    final_integral = sympy.integrate(integral_z, (theta, 0, theta_limit_upper))

    # 4. Display the final equation and answer
    print("Step 3: The result of the outermost integral (with respect to theta) is the final answer.")
    # The numbers in the final equation are the result of the z-integral and the limits for theta.
    print(f"The final calculation is: ∫ from 0 to 2*pi ({integral_z}) d(theta) = ({integral_z}) * (2*pi - 0)")
    print(f"\nThe value of the integral is: {final_integral}")

solve_cone_integral()
<<<108*pi/35>>>
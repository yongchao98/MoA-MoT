import sympy

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2*(x^2+y^2) over a cone.
    The cone has a base radius of 3 and a height of 2.
    """
    # Define the symbols for cylindrical coordinates
    r, theta, z = sympy.symbols('r theta z')

    # The function in cylindrical coordinates is f = z^2 * r^2
    # The integrand is f * r (due to the Jacobian for cylindrical coordinates)
    integrand = z**2 * r**3

    # Define the limits of integration
    # r: from 0 to the line defining the cone's side, r = 3*(1 - z/2)
    r_limits = (r, 0, 3 * (1 - z/2))
    # z: from 0 to the height of the cone, 2
    z_limits = (z, 0, 2)
    # theta: a full circle, 0 to 2*pi
    theta_limits = (theta, 0, 2 * sympy.pi)

    # Perform the integration step-by-step
    # 1. Integrate with respect to r
    integral_r = sympy.integrate(integrand, r_limits)
    # 2. Integrate the result with respect to z
    integral_z = sympy.integrate(integral_r, z_limits)
    # 3. Integrate the result with respect to theta
    final_result = sympy.integrate(integral_z, theta_limits)

    # Extract the numbers from the final symbolic result for clear output
    # The result is expected to be a fraction multiplied by pi
    if sympy.pi in final_result.args:
        coefficient = final_result / sympy.pi
        num, den = sympy.fraction(coefficient)
        
        print("The integral of f(x,y,z) = z^2*(x^2+y^2) over the cone is calculated.")
        print("The integral is set up in cylindrical coordinates as:")
        print("Integral from 0 to 2*pi d(theta) * Integral from 0 to 2 dz * Integral from 0 to 3*(1-z/2) dr of (z^2 * r^3)")
        print("\n--- Final Result ---")
        print(f"The final equation is: Integral = ({num} * pi) / {den}")
        # Output the numbers in the final equation
        print(f"The numbers in the final equation are: numerator = {num}, denominator = {den}.")
        print(f"The symbolic result is: {final_result}")
        print(f"The numeric value is approximately: {final_result.evalf()}")
    else:
        print(f"The result of the integration is: {final_result}")

if __name__ == '__main__':
    solve_cone_integral()
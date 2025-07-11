import sympy

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2*(x^2+y^2) over a cone
    with base radius 3 and height 2, using symbolic integration.
    """
    # Define the symbols for cylindrical coordinates
    r, theta, z = sympy.symbols('r theta z')

    # The function f in cylindrical coordinates is z^2 * r^2.
    # The volume element dV is r * dr * d(theta) * dz.
    # The integrand is f * r.
    integrand = z**2 * r**3

    # Define the integration limits for the cone
    # Height H=2, Radius R=3
    # z from 0 to H
    z_limits = (z, 0, 2)
    # theta for a full circle
    theta_limits = (theta, 0, 2 * sympy.pi)
    # r depends on z: r = R/H * (H - z) = 3/2 * (2-z)
    r_limits = (r, 0, (sympy.S(3)/2) * (2 - z))

    # Explain the setup
    print("This script calculates the integral of f(x,y,z) = z^2(x^2+y^2) over a cone.")
    print("The integral is solved in cylindrical coordinates (r, theta, z).")
    print(f"Integrand: {integrand}")
    print("Integration Limits:")
    print(f"z: from {z_limits[1]} to {z_limits[2]}")
    print(f"theta: from {theta_limits[1]} to {theta_limits[2]}")
    print(f"r: from {r_limits[1]} to {r_limits[2]}")
    print("\nCalculating...")

    # Perform the integration step-by-step
    # 1. Integrate with respect to r
    integral_r = sympy.integrate(integrand, r_limits)

    # 2. Integrate the result with respect to theta
    integral_theta = sympy.integrate(integral_r, theta_limits)

    # 3. Integrate the result with respect to z
    final_integral_value = sympy.integrate(integral_theta, z_limits)

    # Extract the rational part of the answer to display the numbers
    rational_part = final_integral_value / sympy.pi
    num, den = sympy.fraction(rational_part)

    # Print the final result and the numbers in the equation
    print("\n--- Final Result ---")
    print("The final equation for the integral (I) is:")
    print(f"I = ({num} * pi) / {den}")
    
    # Print the numerical value for reference
    numerical_value = final_integral_value.evalf()
    print(f"\nThe numerical value is approximately: {numerical_value}")
    
    return numerical_value

if __name__ == '__main__':
    result = solve_cone_integral()
    # The final output tag, as requested.
    # print(f"\n<<<{result}>>>")

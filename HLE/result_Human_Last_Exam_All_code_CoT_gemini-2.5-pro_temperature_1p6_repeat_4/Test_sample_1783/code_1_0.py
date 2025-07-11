import sympy

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2*(x^2+y^2) over a cone.
    """
    # Define the symbolic variables
    r, theta, z = sympy.symbols('r theta z')
    pi = sympy.pi

    # The function f(x,y,z) = z^2*(x^2+y^2) in cylindrical coordinates is f = z^2*r^2.
    # The integrand includes the Jacobian 'r' from the coordinate transformation.
    integrand = z**2 * r**3

    # Define the limits of integration for the cone
    # Height H=2, Radius R=3
    # z: 0 to 2
    # r: 0 to 3*(1 - z/2)
    # theta: 0 to 2*pi
    limit_r = (r, 0, 3 * (1 - z/2))
    limit_z = (z, 0, 2)
    limit_theta = (theta, 0, 2 * pi)

    # We express the full integral symbolically before calculating
    # The numbers in the equation are:
    # integrand coefficients and powers: 2, 3
    # r limits: 0, 3, 1, 2
    # z limits: 0, 2
    # theta limits: 0, 2*pi
    full_integral = sympy.Integral(integrand, limit_r, limit_z, limit_theta)

    # Calculate the integral
    result = sympy.integrate(integrand, limit_r, limit_z, limit_theta)

    # Print the final equation setup and the result
    print("The integral to solve is:")
    # Using unicode for a cleaner integral representation
    print(f"∫(from {limit_theta[1]} to {limit_theta[2]}) dθ ∫(from {limit_z[1]} to {limit_z[2]}) dz ∫(from {limit_r[1]} to {limit_r[2]}) ({integrand}) dr")
    
    print("\nResult of the integration:")
    print(result)

if __name__ == '__main__':
    solve_cone_integral()
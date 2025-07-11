import sympy
from sympy import integrate, pi, Integral

def solve_cone_integral():
    """
    This function calculates the integral of f(x,y,z) = z^2*(x^2+y^2)
    over the volume of a cone with base radius 3 and height 2.
    """
    # 1. Define the symbols for cylindrical coordinates
    r, theta, z = sympy.symbols('r theta z')

    # 2. Define the integrand in cylindrical coordinates
    # The function f is z^2 * r^2. The volume element dV is r*dr*dtheta*dz.
    # The integrand is f * r.
    integrand = z**2 * r**3

    # 3. Define the limits of integration for the cone
    z_min, z_max = 0, 2
    theta_min, theta_max = 0, 2 * pi
    # The radius r is a function of z: r = 3 - (3/2)*z
    r_min, r_max = 0, 3 - (sympy.S(3)/2) * z

    # 4. Set up the integral expression to show the problem setup
    integral_expression = Integral(integrand, (r, r_min, r_max), (z, z_min, z_max), (theta, theta_min, theta_max))
    print("The integral to be solved in cylindrical coordinates is:")
    print(integral_expression)
    print("-" * 30)

    # 5. Calculate the integral step-by-step
    # First, integrate with respect to r
    integral_r = integrate(integrand, (r, r_min, r_max))
    
    # Second, integrate the result with respect to z
    # This gives the result of the double integral over r and z
    integral_rz = integrate(integral_r, (z, z_min, z_max))
    
    # Finally, integrate with respect to theta. Since the integrand
    # does not depend on theta, this is equivalent to multiplying by 2*pi.
    final_result = integrate(integral_rz, (theta, theta_min, theta_max))

    # 6. Print the final equation and the result, showing the numbers involved
    print("The result is the product of the integral over r and z, and the integral over theta.")
    print(f"The final equation is: ({integral_rz}) * (2*pi) = {final_result}")

if __name__ == '__main__':
    solve_cone_integral()

<<<108*pi/35>>>
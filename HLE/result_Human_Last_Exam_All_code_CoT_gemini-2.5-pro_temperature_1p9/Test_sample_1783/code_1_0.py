import sympy

def solve_cone_integral():
    """
    This function calculates the integral of f(x,y,z) = z^2*(x^2+y^2) over a cone
    with base radius 3 and height 2, using symbolic integration in cylindrical coordinates.
    """
    # Define the symbolic variables
    r, theta, z = sympy.symbols('r theta z')
    
    # Define the function and volume element in cylindrical coordinates.
    # The function f becomes z^2 * r^2.
    # The volume element dV is r*dr*dtheta*dz.
    # The full integrand is (z^2 * r^2) * r
    integrand = z**2 * r**3

    # Define the limits of integration
    # For r (radius): from 0 to the edge of the cone, which is a function of z.
    r_min = 0
    r_max = sympy.Rational(3, 2) * (2 - z)
    
    # For theta (angle): a full circle
    theta_min = 0
    theta_max = 2 * sympy.pi
    
    # For z (height): from base to apex
    z_min = 0
    z_max = 2

    # Step 1: Integrate with respect to r
    integral_r = sympy.integrate(integrand, (r, r_min, r_max))

    # Step 2: Integrate the result with respect to theta
    integral_theta = sympy.integrate(integral_r, (theta, theta_min, theta_max))

    # Step 3: Integrate the final result with respect to z
    final_integral = sympy.integrate(integral_theta, (z, z_min, z_max))

    # The result is a symbolic expression. We extract the numbers to print the final equation.
    if final_integral.is_constant() and sympy.pi in final_integral.args:
        # The result is of the form C * pi, where C is a rational number.
        coefficient = final_integral / sympy.pi
        num = coefficient.p
        den = coefficient.q
        
        print("The final value of the integral is given by the equation:")
        # Output each number in the final equation
        print(f"I = ({num} * pi) / {den}")
        print(f"\nAs a floating point number, this is approximately: {final_integral.evalf()}")
    else:
        # Fallback for unexpected format
        print(f"The final result of the integral is: {final_integral}")

solve_cone_integral()
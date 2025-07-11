import sympy

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2(x^2+y^2) over a cone
    with base radius 3 and height 2, centered at the origin.
    """
    # 1. Define the symbols for cylindrical coordinates
    r, theta, z = sympy.symbols('r theta z')
    pi = sympy.pi

    # 2. Define the integrand in cylindrical coordinates
    # Original function f(x,y,z) = z^2 * (x^2 + y^2)
    # In cylindrical coordinates, f(r,theta,z) = z^2 * r^2
    # The volume element dV is r*dr*dtheta*dz
    # The full integrand is f * r = z^2 * r^3
    integrand = z**2 * r**3

    # 3. Define the limits of integration
    # Height H=2, Radius R=3
    # The radius of the cone at height z is r(z) = R * (1 - z/H) = 3 * (1 - z/2)
    limit_r = (r, 0, 3 * (1 - z/2))
    limit_theta = (theta, 0, 2*pi)
    limit_z = (z, 0, 2)

    # 4. Compute the triple integral by integrating from the inside out
    integral_over_r = sympy.integrate(integrand, limit_r)
    integral_over_theta = sympy.integrate(integral_over_r, limit_theta)
    final_result = sympy.integrate(integral_over_theta, limit_z)

    # 5. Print the result in a structured way
    print("The integral of f(x,y,z) = z^2(x^2+y^2) over the specified cone is calculated as follows:")
    
    # Extract the components of the final fraction for clear output
    p, q = final_result.as_numer_denom()
    
    # The numerator contains pi, so we separate its coefficient
    p_coeff = p / pi
    
    print(f"Final Equation: ({p_coeff} * \u03C0) / {q}")
        
    # Print the numerical value for practical use
    numerical_value = final_result.evalf()
    print(f"Numerical value: {numerical_value}")

if __name__ == '__main__':
    solve_cone_integral()
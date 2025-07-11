import sympy
from sympy import pi, integrate

def solve_cone_integral():
    """
    Calculates the integral of f(x,y,z) = z^2(x^2+y^2) over a cone
    with base radius 3 and height 2.
    """
    # Define the symbols for cylindrical coordinates
    r, theta, z = sympy.symbols('r theta z')

    # Define the function in cylindrical coordinates: f = z^2 * r^2
    # The integrand includes the Jacobian 'r' for the volume element: z^2 * r^3
    integrand = z**2 * r**3

    # Define the geometry of the cone
    R = 3  # Radius of the base
    H = 2  # Height of the cone

    # The upper limit for the radius 'r' depends on 'z'
    r_limit_upper = R * (1 - z / H)

    # Print the description of the integral being solved
    print("The integral to solve is ∫∫∫_V f(x,y,z) dV where f(x,y,z) = z^2(x^2+y^2).")
    print(f"The volume V is a cone with base radius {R} and height {H}.")
    print("\nIn cylindrical coordinates, the integral is set up as:")
    print(f"I = ∫(θ=0 to 2π) dθ ⋅ ∫(z=0 to {H}) dz ⋅ ∫(r=0 to {R}(1-z/{H})) [ {integrand} ] dr")
    
    # Perform the integration symbolically
    # 1. Integrate with respect to r
    integral_r = integrate(integrand, (r, 0, r_limit_upper))
    
    # 2. Integrate the result with respect to z
    integral_z = integrate(integral_r, (z, 0, H))
    
    # 3. Integrate the result with respect to theta
    final_result = integrate(integral_z, (theta, 0, 2*pi))

    # Deconstruct the final answer to show its components as requested
    if final_result.is_rational:
        num_val = final_result
        den_val = 1
        pi_str = ""
    else:
        num_val, den_val = (final_result / pi).as_numer_denom()
        pi_str = " * π"

    print("\nAfter calculation, the final result is:")
    print(f"Final Equation: Integral = ({num_val}{pi_str}) / {den_val}")
    print(f"Which corresponds to the numerical value: {final_result.evalf()}")

solve_cone_integral()
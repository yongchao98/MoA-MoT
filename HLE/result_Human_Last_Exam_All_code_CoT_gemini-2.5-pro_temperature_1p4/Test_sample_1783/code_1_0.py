import sympy as sp

def solve_integral():
    """
    Calculates the integral of f(x,y,z) = z^2(x^2+y^2) over a cone.
    The cone has a base of radius 3 at z=0 and a vertex at height z=2.
    The calculation is performed step-by-step using cylindrical coordinates.
    """
    # Define symbolic variables
    r, theta, z = sp.symbols('r theta z')
    pi = sp.pi

    # Define the function in cylindrical coordinates and the volume element (Jacobian)
    f = z**2 * r**2
    jacobian = r
    integrand = f * jacobian

    # Define the limits of integration
    # For r: 0 to (3/2)*(2-z)
    # For z: 0 to 2
    # For theta: 0 to 2*pi
    r_upper_limit = sp.Rational(3, 2) * (2 - z)
    r_limits = (r, 0, r_upper_limit)
    z_limits = (z, 0, 2)
    theta_limits = (theta, 0, 2 * pi)

    # --- Print the setup ---
    print("The task is to compute the integral of f(x,y,z) = z^2*(x^2+y^2) inside a cone.")
    print("We switch to cylindrical coordinates.")
    print(f"The function becomes: f(r, z) = {f}")
    print(f"The volume element dV becomes: {jacobian} dr d(theta) dz")
    print(f"The integrand is: {integrand}")
    print(f"The limits are:")
    print(f"  r: from 0 to {r_upper_limit}")
    print(f"  z: from 0 to 2")
    print(f"  theta: from 0 to 2*pi\n")

    # --- Step-by-step integration ---
    # 1. Integrate with respect to r
    print("Step 1: Integrating with respect to r...")
    integral_r = sp.integrate(integrand, r_limits)
    print(f"Result after r-integration: {integral_r}\n")

    # 2. Integrate with respect to z
    print("Step 2: Integrating the result with respect to z...")
    integral_z = sp.integrate(integral_r, z_limits)
    print(f"Result after z-integration: {integral_z}\n")

    # 3. Integrate with respect to theta
    print("Step 3: Integrating the result with respect to theta...")
    final_result = sp.integrate(integral_z, theta_limits)
    print(f"Result after theta-integration: {final_result}\n")
    
    # --- Final Answer Output ---
    # Extract the numerical components of the final symbolic answer
    if isinstance(final_result, sp.Mul) and pi in final_result.args:
        # The result is in the form of C * pi
        coeff = final_result.as_coeff_Mul(pi)[0]
        if isinstance(coeff, sp.Rational):
            num_val = coeff.p
            den_val = coeff.q
            pi_symbol = 'pi'
            
            print("The final result is an equation in the form (numerator * pi) / denominator.")
            print("Final Equation:")
            print(f"({num_val} * {pi_symbol}) / {den_val}")
    else:
        # If the result is not a multiple of pi, just print it
        print(f"Final symbolic answer: {final_result}")

    # Print the numerical value
    numerical_value = final_result.evalf()
    print(f"\nNumerical value: {numerical_value}")

solve_integral()
<<<9.8>>>
import sympy as sp

def evaluate_1s_integral():
    """
    This function evaluates the integral <phi_1s | 1/r | phi_1s> for a 
    1s Slater-type orbital using symbolic mathematics with the sympy library.
    It demonstrates the process step-by-step, from normalization to the final result.
    """
    # Define the necessary symbolic variables for spherical coordinates (r, theta, phi)
    # and the orbital exponent (zeta).
    r, theta, phi = sp.symbols('r theta phi', real=True, positive=True)
    zeta = sp.symbols('zeta', real=True, positive=True)

    # Use the infinity symbol from sympy for integration limits.
    oo = sp.oo

    # Step 1: Define the unnormalized 1s Slater-type orbital (STO).
    # The functional form is an exponential decay.
    phi_1s_unnormalized = sp.exp(-zeta * r)

    print("Step 1: The unnormalized 1s Slater-type orbital (STO) is:")
    print(f"    \u03C6_unnormalized(r) = {phi_1s_unnormalized}\n")

    # Step 2: Calculate the normalization constant, N.
    # The normalization condition is Integral(|N * phi|^2 dV) = 1 over all space.
    # The volume element in spherical coordinates is dV = r^2 * sin(theta) * dr * dtheta * dphi.
    integrand_for_norm = phi_1s_unnormalized**2 * r**2 * sp.sin(theta)
    
    # We solve the integral of |phi_unnormalized|^2 dV.
    norm_integral_val = sp.integrate(integrand_for_norm, (phi, 0, 2*sp.pi), (theta, 0, sp.pi), (r, 0, oo))

    # N is the reciprocal of the square root of this integral.
    N = sp.sqrt(1 / norm_integral_val)
    
    print("Step 2: Calculate the normalization constant N.")
    print(f"    The integral of |\u03C6_unnormalized|^2 over all space evaluates to: {norm_integral_val}")
    print(f"    The normalization constant N is therefore: {N}\n")

    # Step 3: Define the normalized 1s STO.
    phi_1s_normalized = N * phi_1s_unnormalized
    print("Step 3: The normalized 1s STO is:")
    print(f"    \u03C6_1s(r) = {phi_1s_normalized}\n")

    # Step 4: Evaluate the integral <phi_1s | 1/r | phi_1s>.
    # The operator is 1/r.
    operator = 1 / r
    
    # The integrand is |phi_1s|^2 * (1/r) * dV.
    final_integrand = phi_1s_normalized**2 * operator * r**2 * sp.sin(theta)

    print("Step 4: Set up the final integral for <\u03C6_1s| 1/r |\u03C6_1s>.")
    print(f"    The integrand is |{phi_1s_normalized}|^2 * ({operator}) * r^2 * sin(\u03B8)\n")

    # Perform the final integration over all space.
    final_result = sp.integrate(final_integrand, (phi, 0, 2*sp.pi), (theta, 0, sp.pi), (r, 0, oo))

    # Step 5: Display the final result.
    # The final equation is <phi_1s | (1/r) | phi_1s> = result.
    op_num = 1 # The number in the operator term 1/r
    
    print("Step 5: The result of the integration is:")
    print(f"    <\u03C6_1s| ({op_num}/r) |\u03C6_1s> = {final_result}")
    print("\nConclusion: The expectation value of the 1/r operator for a 1s Slater orbital is simply the orbital exponent, zeta.")
    print("For example, for a hydrogen atom, where \u03B6 = 1, the value is 1.")

if __name__ == '__main__':
    evaluate_1s_integral()
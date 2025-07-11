import sympy

def evaluate_slater_integral():
    """
    Symbolically evaluates the integral <phi_1s | 1/r | phi_1s> for a Slater-Type Orbital (STO)
    using the sympy library.
    """
    
    # --- Step 1: Define symbolic variables ---
    # r: radial distance
    # theta: polar angle
    # phi: azimuthal angle
    # zeta: the orbital exponent, a positive real number
    r, theta, phi, zeta = sympy.symbols('r theta phi zeta', real=True, positive=True)

    # --- Step 2: Define the 1s Slater Orbital ---
    # The unnormalized 1s Slater orbital has the functional form exp(-zeta*r).
    # The full orbital is phi(r) = N * exp(-zeta*r), where N is the normalization constant.
    phi_1s_unnormalized = sympy.exp(-zeta * r)
    
    print("This script evaluates the integral <phi_1s | 1/r | phi_1s>.")
    print(f"The unnormalized 1s Slater orbital is defined as: {phi_1s_unnormalized}\n")
    
    # The volume element in spherical coordinates is r^2 * sin(theta)
    volume_element = r**2 * sympy.sin(theta)

    # --- Step 3: Calculate the Normalization Constant N ---
    # The normalization condition is that the integral of |phi|^2 over all space is 1.
    # Integral( (N * exp(-zeta*r))^2 * dV ) = 1
    print("--- Calculating the Normalization Constant N ---")
    integrand_norm = phi_1s_unnormalized**2 * volume_element
    
    # We integrate over all space: r from 0 to infinity, theta from 0 to pi, phi from 0 to 2*pi
    norm_integral_val = sympy.integrate(integrand_norm, (phi, 0, 2*sympy.pi), (theta, 0, sympy.pi), (r, 0, sympy.oo))
    
    print(f"The normalization integral <phi_unnormalized|phi_unnormalized> is: {norm_integral_val}")
    
    # From N^2 * norm_integral_val = 1, we find N.
    N_squared = 1 / norm_integral_val
    N = sympy.sqrt(N_squared)
    phi_1s_normalized = N * phi_1s_unnormalized
    
    print(f"The normalization constant squared, N^2, is: {N_squared}")
    print(f"The normalized 1s Slater orbital is: {phi_1s_normalized}\n")
    
    # --- Step 4: Evaluate the expectation value <phi_1s | 1/r | phi_1s> ---
    print("--- Evaluating the Integral <phi_1s | 1/r | phi_1s> ---")
    
    # The operator is 1/r
    operator = 1 / r
    
    # The integrand for the expectation value is (phi_normalized*) * operator * (phi_normalized) * dV
    integrand_expectation = phi_1s_normalized.conjugate() * operator * phi_1s_normalized * volume_element
    
    print("The integrand is: (phi_normalized*) * (1/r) * (phi_normalized) * r^2 * sin(theta)")
    
    # Perform the final triple integration
    expectation_value = sympy.integrate(integrand_expectation, (phi, 0, 2*sympy.pi), (theta, 0, sympy.pi), (r, 0, sympy.oo))
    
    # --- Step 5: Display the Final Result ---
    print("\n--- Final Result ---")
    final_equation_str = f"< phi_1s | 1/r | phi_1s > = {expectation_value}"
    print(final_equation_str)

if __name__ == '__main__':
    evaluate_slater_integral()
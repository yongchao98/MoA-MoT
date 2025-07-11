import math

def approximate_integral():
    """
    This function calculates the analytical approximation for the integral
    I(epsilon) for small epsilon, based on the method described above.
    """
    
    # Parameters from the leading term p(x) ≈ c * x^alpha near x=0.
    # For p(x) = 9.0*x^5.0 + 5.0*x^6.0 + 9.0*x^8.0, we have:
    c = 9.0
    alpha = 5.0
    
    # The exponent for epsilon in the approximation formula is p = -(alpha-1)/alpha.
    exponent = -(alpha - 1.0) / alpha
    
    # The coefficient in the approximation is C = (pi / (alpha * sin(pi/alpha) * c^(1/alpha))).
    pi_over_alpha = math.pi / alpha
    sin_pi_over_alpha = math.sin(pi_over_alpha)
    c_to_one_over_alpha = c**(1.0/alpha)
    
    coefficient = (math.pi / alpha) / (sin_pi_over_alpha * c_to_one_over_alpha)
    
    # Print the final derived formula with the calculated numerical values.
    # The final equation includes all the necessary numbers.
    print("The analytical approximation for I(\u03B5) in the small \u03B5 regime is:")
    print(f"I(\u03B5) \u2248 C \u00B7 \u03B5\u1D56")
    print(f"where the exponent is p = {exponent}")
    print(f"and the coefficient is C = {coefficient}")
    print("\nThus, the final formula is:")
    print(f"I(\u03B5) \u2248 {coefficient} \u00B7 \u03B5^({exponent})")

# Execute the function to display the result.
approximate_integral()
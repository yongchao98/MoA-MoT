import math

def estimate_alpha_from_theory(n_values):
    """
    This function explains the theoretical derivation of alpha based on
    the relationship between a polynomial's degree and its maximum derivative.
    """
    print("--- Theoretical Derivation of alpha ---")
    print("The minimum degree d_n of the polynomial is determined by how fast it must change.")
    print("The polynomial p_n(x) must transition from <= 1 at x=n^2 to >= 2 at x=n^2+1.")

    print("\nStep 1: Scale the coordinates.")
    print("We map the interval [n^2+1, n^10] to [-1, 1].")
    print("The point x=n^2 is mapped to v(n^2) approx -1 - 2/n^8.")
    print("The point x=n^2+1 is mapped to v(n^2+1) = -1.")
    
    print("\nStep 2: Use the Mean Value Theorem.")
    print("Let P(v) be the polynomial in the scaled coordinate v.")
    print("P(-1) >= 2 and P(-1 - 2/n^8) <= 1.")
    print("This implies the derivative P'(c) must be at least (2-1) / (2/n^8) = n^8 / 2.")

    print("\nStep 3: Use Markov's Inequality for polynomials.")
    print("This inequality bounds the derivative: |P'(v)| <= M * d_n^2, where M is the max value on [-1,1].")
    print("Here, M is approximately 3 (from the condition p_n(i) in [2,3]).")

    print("\nStep 4: Combine the bounds to find the relation for d_n.")
    print("d_n^2 * M >= n^8 / 2  =>  d_n >= n^4 / sqrt(2*M)")
    print("This shows that d_n grows proportionally to n^4. Thus, alpha = 4.")

    print("\n--- Numerical Example ---")
    for n in n_values:
        # The gap in the scaled coordinate system
        delta = 2 / (n**8)
        
        # Lower bound on the derivative from Mean Value Theorem
        min_derivative = 1 / delta
        
        # Maximum value of the polynomial on the interval [-1, 1]
        M = 3.0
        
        # Lower bound on the degree squared from Markov's Inequality
        d_n_squared_lower_bound = min_derivative / M
        
        # Lower bound on the degree
        d_n_lower_bound = math.sqrt(d_n_squared_lower_bound)
        
        print(f"\nFor n = {n}:")
        print(f"The required change must happen over a scaled distance delta = 2/n^8 = {delta:.2e}")
        final_equation = f"d_{n} >= sqrt( (n^8 / 2) / M ) = sqrt( ({n**8 / 2:.2e}) / {M} ) = {d_n_lower_bound:.2e}"
        print(f"The final equation for the lower bound on the degree is:")
        print(final_equation)
        
        # The asymptotic form is d_n = C * n^alpha
        # d_n_lower_bound = (1/sqrt(2*M)) * n^4
        # So we can see the n^4 dependence.
        
# Run the explanation and calculation for a sample value of n
estimate_alpha_from_theory([10])

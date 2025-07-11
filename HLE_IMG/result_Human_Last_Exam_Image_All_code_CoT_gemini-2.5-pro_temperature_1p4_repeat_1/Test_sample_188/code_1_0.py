import sympy

def compute_chromatic_polynomial():
    """
    Computes the chromatic polynomial for the given graph and prints the result.
    The derivation is based on the structural decomposition of the graph.
    """
    # Define k as a symbolic variable for the polynomial
    k = sympy.Symbol('k')

    # Based on the step-by-step derivation, the chromatic polynomial in factored form is:
    # P(G, k) = k * (k-1) * (k-2) * (k^2 - 4*k + 5)
    poly_factored = k * (k - 1) * (k - 2) * (k**2 - 4 * k + 5)

    # Use the sympy library to expand the polynomial into its standard form
    poly_expanded = sympy.expand(poly_factored)

    # Extract the integer coefficients from the expanded polynomial object
    coeffs = [int(c) for c in sympy.Poly(poly_expanded, k).all_coeffs()]
    c5, c4, c3, c2, c1 = coeffs

    # Print the final equation, explicitly showing each coefficient as requested.
    # The format is chosen for readability.
    print("The chromatic polynomial of the graph is:")
    print(f"P(G, k) = {c5}*k^5 - {abs(c4)}*k^4 + {c3}*k^3 - {abs(c2)}*k^2 + {c1}*k")

compute_chromatic_polynomial()
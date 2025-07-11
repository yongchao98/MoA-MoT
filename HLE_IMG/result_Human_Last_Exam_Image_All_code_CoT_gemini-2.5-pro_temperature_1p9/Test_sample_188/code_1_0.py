import sympy

def solve_chromatic_polynomial():
    """
    Computes and displays the chromatic polynomial for the given graph
    using the deletion-contraction method.
    """
    # Define k as a symbolic variable for our polynomial
    k = sympy.Symbol('k')

    # Step 1: Define the polynomials for G-e and G/e based on our derivation.
    # Edge e is (1,5).
    # P(G-e, k) = k(k-1)^2(k-2)^2
    p_g_minus_e = k * (k - 1)**2 * (k - 2)**2

    # P(G/e, k) = P(K4, k) = k(k-1)(k-2)(k-3)
    p_g_contract_e = k * (k - 1) * (k - 2) * (k - 3)

    # Step 2: Apply the deletion-contraction formula: P(G, k) = P(G-e, k) - P(G/e, k)
    chromatic_poly = p_g_minus_e - p_g_contract_e

    # Step 3: Expand the resulting polynomial to its standard form.
    expanded_poly = sympy.expand(chromatic_poly)
    
    # Step 4: Format the polynomial for clear output, showing all terms and coefficients.
    # Sympy's string representation uses Python's syntax (e.g., ** for power).
    # We will replace it with conventional mathematical notation (^ for power).
    poly_str = str(expanded_poly).replace('**', '^').replace('*', '')

    print("The chromatic polynomial P(k) is calculated using the deletion-contraction method.")
    print("Formula: P(G, k) = P(G-e, k) - P(G/e, k) for edge e=(1,5)")
    print(f"P(G-e, k) = k(k-1)^2(k-2)^2")
    print(f"P(G/e, k) = k(k-1)(k-2)(k-3)")
    print("\nResulting Chromatic Polynomial:")
    # The derived factored form: k * (k-1) * (k-2) * (k^2 - 4k + 5)
    print("Factored Form: P(k) = k * (k - 1) * (k - 2) * (k^2 - 4k + 5)")
    print("\nExpanded Form:")
    # The numbers in the final equation are 1, -7, 19, -23, 10
    final_equation = f"P(k) = {poly_str}"
    print(final_equation)


solve_chromatic_polynomial()
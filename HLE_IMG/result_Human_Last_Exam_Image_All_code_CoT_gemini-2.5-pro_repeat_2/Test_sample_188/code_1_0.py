import sympy

def solve():
    """
    Computes the chromatic polynomial for the given graph using the deletion-contraction method.
    """
    k = sympy.Symbol('k')

    # Step 1: Define the chromatic polynomials of the component graphs derived from deletion-contraction.

    # P(G/e, k): Chromatic polynomial of G after contracting edge e=(4,5).
    # The resulting graph is a K4-e (a complete graph on 4 vertices minus one edge).
    # Its polynomial is k(k-1)(k-2)^2.
    P_G_contract_e = k * (k - 1) * (k - 2)**2

    # P(G-e, k): Chromatic polynomial of G after deleting edge e=(4,5).
    # This graph G' requires another deletion-contraction step on edge e'=(3,4).
    # P(G', k) = P(G'-e', k) - P(G'/e', k)

    # P(G'-e', k): G-e-e' is a C4 (1-2-3-5) with a pendant edge (2,4).
    # Its polynomial is k * (k-1)^2 * (k^2 - 3*k + 3).
    P_G_prime_delete_e_prime = k * (k - 1)**2 * (k**2 - 3*k + 3)

    # P(G'/e', k): G-e contracted by e'=(3,4) is a C4.
    # The polynomial for C4 is k^4 - 4k^3 + 6k^2 - 3k, or k(k-1)(k^2-3k+3)
    P_G_prime_contract_e_prime = k * (k - 1) * (k**2 - 3*k + 3)
    
    # Calculate P(G-e, k)
    P_G_delete_e = P_G_prime_delete_e_prime - P_G_prime_contract_e_prime

    # Step 2: Calculate the final chromatic polynomial for the original graph G.
    # P(G, k) = P(G-e, k) - P(G/e, k)
    P_G = P_G_delete_e - P_G_contract_e

    # Step 3: Expand the polynomial to get the final expression.
    final_poly = sympy.expand(P_G)
    
    # Step 4: Format and print the output.
    print("The chromatic polynomial P(k) for the graph is:")
    
    # Extract coefficients and print in the specified format
    poly_coeffs = final_poly.as_poly(k).all_coeffs()
    degree = len(poly_coeffs) - 1
    
    equation_parts = []
    for i, coeff in enumerate(poly_coeffs):
        power = degree - i
        if coeff == 0:
            continue
        
        # Format coefficient
        if power > 0 and abs(coeff) == 1:
            coeff_str = "" if coeff == 1 else "-"
        else:
            coeff_str = str(coeff)
        
        # Format variable part
        if power == 0:
            var_str = ""
        elif power == 1:
            var_str = "k"
        else:
            var_str = f"k**{power}"
        
        # Combine coefficient and variable
        if var_str and coeff_str and coeff_str != "-":
            term = f"{coeff_str}*{var_str}"
        else:
            term = f"{coeff_str}{var_str}"

        equation_parts.append(term)
        
    final_equation = " + ".join(equation_parts).replace("+ -", "- ")
    
    print(f"P(k) = {final_poly}")
    print("\nThe final equation with each coefficient is:")
    print(f"P(k) = {final_equation}")

solve()
import sympy

def solve():
    """
    Computes the chromatic polynomial for the given graph using the deletion-contraction method.
    """
    k = sympy.Symbol('k')

    # Step 1: Define the chromatic polynomial for the graph G' = G - e, where e = (3,4).
    # As derived in the plan, P(G', k) = k * (k-1) * (k-2)**3
    P_G_minus_e = k * (k - 1) * (k - 2)**3

    # Step 2: Define the chromatic polynomial for the graph G'' = G / e.
    # G'' is a 4-cycle, P(C4, k) = (k-1)**4 + (k-1)
    # Let's expand it: k**4 - 4*k**3 + 6*k**2 - 4*k + 1 + k - 1 = k**4 - 4*k**3 + 6*k**2 - 3*k
    P_G_contract_e = k**4 - 4*k**3 + 6*k**2 - 3*k

    # Step 3: Apply the deletion-contraction formula P(G, k) = P(G', k) - P(G'', k)
    P_G = sympy.expand(P_G_minus_e - P_G_contract_e)

    # Step 4: Format the output string to show each term of the polynomial.
    # The result is P_G = k**5 - 7*k**4 + 19*k**3 - 23*k**2 + 10*k
    
    # Let's build the output string from the sympy result
    poly_terms = []
    # as_ordered_terms() gives terms from highest power to lowest
    for power, coeff in P_G.as_poly(k).all_terms():
        # Format coefficient
        if coeff == 1 and power != 0:
            coeff_str = ""
        elif coeff == -1 and power != 0:
            coeff_str = "-"
        else:
            coeff_str = str(coeff)

        # Format power
        if power == 0:
            var_str = ""
        elif power == 1:
            var_str = "k"
        else:
            var_str = f"k^{power}"
        
        if coeff != 0:
            if coeff > 0:
                sign = " + "
            else:
                sign = " - "
                # remove negative sign from coeff_str if it exists
                if coeff_str.startswith('-'):
                    coeff_str = coeff_str[1:]

            term_str = f"{coeff_str}{var_str}"
            poly_terms.append((sign, term_str))

    # Join the terms into the final equation string
    # First term doesn't need a leading sign if it's positive
    if poly_terms[0][0] == " + ":
        final_string = "P(G, k) = " + poly_terms[0][1]
    else: # if it's negative
        final_string = "P(G, k) = -" + poly_terms[0][1]
    
    for sign, term in poly_terms[1:]:
        final_string += sign + term

    print("The chromatic polynomial P(G, k) is calculated using deletion-contraction.")
    print("P(G, k) = P(G - e, k) - P(G / e, k)")
    print(f"P(G - e, k) = k(k-1)(k-2)^3 = {sympy.expand(P_G_minus_e)}")
    print(f"P(G / e, k) = k^4 - 4k^3 + 6k^2 - 3k")
    print("\nThe final chromatic polynomial is:")
    print(final_string)

solve()
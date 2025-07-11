import sympy

def solve_extremal_problem():
    """
    Solves the theoretical graph theory problem and prints the answer in the required format.
    """

    # Part (a): True or false: If G is not a union of K_2's, then ex(n; G, K_{1,t}-ind) = Θ(n).
    # Analysis: This is False. A counterexample can be constructed. For G=P_4 (which is not a union of K_2's),
    # consider a complete split graph H with a clique C of size k=n-t+1 and an independent set I of size n-k=t-1,
    # with all edges between C and I. This graph is P_4-free (as it's a cograph) and K_{1,t}-induced-free,
    # but has O(n^2) edges.
    answer_a = "False"

    # Part (b): If G ~ sK_2 for s >= 2, is ex(n; sK_2, K_{1,t}-ind) independent of n?
    # Analysis: This is True. A graph H with no sK_2 has a maximum matching M of size m < s. The vertices V(M)
    # have size 2m <= 2(s-1). The remaining vertices I = V \ V(M) form an independent set. All edges are
    # incident to V(M). The number of edges within V(M) is at most binom(2(s-1), 2). The number of edges
    # between V(M) and I is at most 2(s-1)(t-1) due to the K_{1,t}-induced-free condition.
    # The total number of edges is bounded by a constant depending only on s and t.
    answer_b = "True"

    # Part (c): For G ~ sK_2, express the upper bound for ex(n; sK_2, K_{1,t}-ind).
    # Analysis: From the logic for (b), the number of edges e(H) is at most
    # f(m) = binom(2m, 2) + 2m(t-1) for m <= s-1. This function increases with m, so we maximize
    # it at m = s-1.
    # Bound = binom(2(s-1), 2) + 2(s-1)(t-1)
    #       = (2s-2)(2s-3)/2 + 2(s-1)(t-1)
    #       = (s-1)(2s-3) + 2(s-1)(t-1)
    #       = (s-1)(2s - 3 + 2t - 2)
    #       = (s-1)(2s + 2t - 5)
    
    s, t = sympy.symbols('s t')
    # Expression for the upper bound
    bound_expr = (s - 1) * (2*s + 2*t - 5)
    # Format for printing, including explicit numbers
    answer_c = "(s - 1)*(2*s + 2*t - 5)"


    # Final output string
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer)

solve_extremal_problem()
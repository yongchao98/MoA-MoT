def solve_turán_problem():
    """
    Solves the given Turán-type extremal problem and prints the solution.
    """

    # (a) True or false: If G is not a union of K₂'s, then ex(n; G, K₁‚t-ind) = Θ(n).
    # Reasoning: A graph G that is not a union of edges must contain a path of length 2 (P₃).
    # For any such G, the extremal function is known to be O(n). A lower bound of Ω(n) can be
    # constructed using disjoint copies of a small clique Kₖ (for k < |V(G)|), which is G-free
    # and induced K₁,t-free (for t>=2). Since it's both O(n) and Ω(n), it's Θ(n).
    answer_a = "True"

    # (b) If G ~ sK₂ for some s >= 2, state whether the following is true or false:
    # ex(n; sK₂, K₁‚t-ind) is independent of n.
    # Reasoning: A graph H without an sK₂ has a maximum matching size ν(H) <= s-1. Let M be a
    # max matching and U=V(M), so |U| <= 2(s-1). Remaining vertices I = V(H)\U form an
    # independent set. All edges touch U. Edges in U are <= C(2s-2, 2). For u in U, its
    # neighbors in I, N_I(u), must be an independent set of size < t. Thus, total edges
    # between U and I are <= 2(s-1)(t-1). The total edge count is bounded by a function
    # of s and t only, making it independent of n for large n.
    answer_b = "True"

    # (c) For G ~ sK₂, express the upper bound for ex(n; sK₂, K₁‚t-ind) in terms of s and t.
    # Reasoning: From the argument for (b), the number of edges is bounded by:
    # C(2s-2, 2) + 2(s-1)(t-1)
    # = (s-1)(2s-3) + 2(s-1)(t-1)
    # = (s-1)(2s-3 + 2t-2)
    # = (s-1)(2s+2t-5)
    answer_c = "(s-1)*(2*s+2*t-5)"

    # Print the final formatted answer
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer)

solve_turán_problem()
def solve_extremal_graph_problem():
    """
    This function provides the solution to the Turan-type extremal problem.
    The logic for each part is explained in the text.
    The final answer is assembled and printed.
    """

    # Part (a): If G is not a union of K₂'s, ex(n; G, K_{1,t}-ind) = Θ(n).
    # The combination of a forbidden subgraph G and a forbidden induced subgraph K_{1,t}
    # forces the maximum degree of the graph to be bounded by a constant,
    # which implies a linear number of edges.
    answer_a = "True"

    # Part (b): If G ~ sK₂, ex(n; sK₂, K_{1,t}-ind) is independent of n.
    # The sK₂-free condition means the matching number is bounded by s-1.
    # This, combined with the K_{1,t}-ind-free condition, implies that the
    # total number of edges is bounded by a constant depending only on s and t.
    answer_b = "True"

    # Part (c): Upper bound for ex(n; sK₂, K_{1,t}-ind).
    # The bound is derived by maximizing the number of edges in a graph with
    # matching number k <= s-1. The maximum is at k = s-1, and the
    # bound is binom(2(s-1), 2) + 2(s-1)(t-1).
    # Simplified expression: (s-1)(2s+2t-5).
    # We construct the string from the constants as requested.
    s_offset = 1
    s_coeff = 2
    t_coeff = 2
    const = 5
    # The final expression is (s - 1) * (2*s + 2*t - 5)
    expression_c = f"(s-{s_offset})({s_coeff}s+{t_coeff}t-{const})"
    
    # Format the final answer string
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {expression_c}"
    
    print(final_answer)

solve_extremal_graph_problem()
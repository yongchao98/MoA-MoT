def solve():
    """
    This function provides the solution to the user's question.
    """

    # Analysis for part (a)
    # The statement is True. A graph that is not a union of K2s contains a component
    # with at least 3 vertices. An argument based on minimum degree shows that if the
    # number of edges is super-linear, the graph must contain the forbidden subgraph G.
    # A construction with disjoint cliques shows the number of edges can be linear.
    # Hence, the function is Theta(n).
    answer_a = "True"

    # Analysis for part (b)
    # The statement is True. A graph with no sK2 has a maximum matching of size at most s-1.
    # Using a structural decomposition based on a maximum matching (V = V_M U I),
    # the number of vertices in V_M is bounded by 2(s-1). The induced K_{1,t}-free property
    # bounds the number of edges connected to the independent set I. This results in a total
    # edge count that depends only on s and t, not n.
    answer_b = "True"

    # Analysis for part (c)
    # Following the argument for (b), we can write down an explicit upper bound.
    # Let k be the size of the max matching, k <= s-1.
    # The number of edges is at most binom(2k, 2) + 2k(t-1).
    # This function is increasing in k, so we maximize it at k = s-1.
    # Bound = binom(2(s-1), 2) + 2(s-1)(t-1)
    #       = (s-1)(2s-3) + 2(s-1)(t-1)
    #       = (s-1)(2s-3 + 2t-2)
    #       = (s-1)(2s+2t-5)
    answer_c_expression = "(s-1)*(2*s+2*t-5)"
    
    # Printing the final answer in the required format
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c_expression}"
    print(final_answer)

solve()
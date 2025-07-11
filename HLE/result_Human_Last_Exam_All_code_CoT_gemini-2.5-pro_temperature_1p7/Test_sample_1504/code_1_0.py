def solve_extremal_graph_problem():
    """
    Solves the three-part graph theory question and prints the answer in the specified format.
    """

    # Part (a): G is not a union of K2's.
    # Reasoning: The function is Omega(n) by constructing a graph of disjoint cliques smaller than G.
    # It is O(n) because the K_{1,t}-induced-free property implies that neighborhoods have bounded
    # independence number. By Ramsey's theorem, this implies the maximum degree of the graph is bounded
    # by a constant, leading to a linear number of edges. Therefore, it's Theta(n).
    answer_a = "True"

    # Part (b): G ~ sK2. Is ex(n; G, K_{1,t}-ind) independent of n?
    # Reasoning: A graph H without an sK_2 matching has a maximum matching size of at most s-1.
    # The vertices V can be partitioned into V(M) (covered by the matching) and I (an independent set).
    # |V(M)| <= 2(s-1). The K_{1,t}-induced-free property restricts the number of edges from any
    # vertex in V(M) to the independent set I. This leads to a total edge count bounded by a
    # constant that only depends on s and t.
    answer_b = "True"

    # Part (c): Upper bound for the sK2 case.
    # Reasoning: Based on the analysis for (b), the number of edges is bounded by the sum of
    # edges within V(M) and edges between V(M) and I.
    # Max edges = C(|V(M)|, 2) + |V(M)| * (t-1), where |V(M)| <= 2(s-1).
    # This gives the expression C(2s-2, 2) + (2s-2)(t-1).
    # Simplified: (s-1)(2s+2t-5). We use the combination form for clarity.
    answer_c = "C(2s-2, 2) + (2s-2)(t-1)"

    # Format the final output string.
    final_answer_string = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    
    print(final_answer_string)

solve_extremal_graph_problem()
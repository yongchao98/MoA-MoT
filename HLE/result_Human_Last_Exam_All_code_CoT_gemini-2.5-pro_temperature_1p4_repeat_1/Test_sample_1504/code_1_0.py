def solve_extremal_graph_problem():
    """
    Solves the Turan-type extremal problem parts (a), (b), and (c).
    The reasoning for each answer is provided in the comments.
    """

    # (a) True or false: If G is not a union of K_2's, then ex(n; G, K_{1,t}-ind) = Theta(n).
    #
    # Reasoning:
    # Let G_0 be a graph that is not a union of K_2's. This means G_0 has a connected component with at least 3 vertices.
    # Lower bound (Ω(n)): Consider a perfect matching on n vertices (if n is even), which is floor(n/2) disjoint K_2's.
    # This graph has floor(n/2) edges, so the number of edges is Omega(n).
    # It is G_0-free because its largest component has 2 vertices, while G_0 has a component with at least 3 vertices.
    # It is K_{1,t}-ind-free for t>=2 because its maximum degree is 1, but an induced K_{1,t} requires a vertex of degree at least t.
    # Upper bound (O(n)): A known result in extremal graph theory states that K_{1,t}-ind-free graphs are structurally close to a disjoint union of cliques.
    # For such a graph to be G_0-free, the size of these cliques must be bounded by the number of vertices in the largest component of G_0.
    # A graph composed of a disjoint union of cliques of bounded size has a total number of edges that is linear in n.
    # Since the function is both Ω(n) and O(n), it is Θ(n).
    answer_a = "True"

    # (b) If G ~ sK_2 for some s >= 2, state whether the following is true or false:
    # ex(n; sK_2, K_{1,t}-ind) is independent of n.
    #
    # Reasoning: "Independent of n" is interpreted as being O(1), i.e., bounded by a constant.
    # Let G' be a graph that is sK_2-free and K_{1,t}-ind-free.
    # sK_2-free means the matching number ν(G') < s, so ν(G') <= s-1.
    # Let M be a maximum matching in G'. Let V(M) be the set of vertices covered by M. |V(M)| <= 2(s-1).
    # The set I = V(G') \ V(M) must be an independent set.
    # Any edge in G' must be incident to at least one vertex in V(M).
    # The number of edges within V(M) is at most C(2(s-1), 2), which is a constant dependent on s.
    # For any vertex v in V(M), its neighbors in I, N(v) ∩ I, form an independent set. Since G' is K_{1,t}-ind-free,
    # the independence number of the subgraph induced by N(v) is less than t. This implies |N(v) ∩ I| <= t-1.
    # The total number of edges between V(M) and I is at most |V(M)| * (t-1) <= 2(s-1)(t-1), a constant.
    # The total number of edges is bounded by C(2s-2, 2) + 2(s-1)(t-1), which does not depend on n.
    answer_b = "True"

    # (c) For G ~ sK_2, express the upper bound for ex(n; sK_2, K_{1,t}-ind) in terms of s and t.
    #
    # Reasoning: From the argument in (b), the upper bound is the sum of the maximum possible edges
    # within V(M) and between V(M) and I.
    # Max edges = C(|V(M)|, 2) + |V(M)|(t-1).
    # This function increases with |V(M)|, so we take the maximum possible size, |V(M)| = 2(s-1).
    # Upper bound = C(2(s-1), 2) + 2(s-1)(t-1)
    #             = (2(s-1))(2(s-1)-1)/2 + 2(s-1)(t-1)
    #             = (s-1)(2s-3) + 2(s-1)(t-1)
    #             = (s-1) * ( (2s-3) + (2t-2) )
    #             = (s-1)(2s+2t-5)
    s = 's'
    t = 't'
    one = 1
    two = 2
    five = 5
    # The required format is an expression. We output the simplified formula as a string.
    # The instruction "output each number in the final equation" is interpreted by showing
    # the final simple integers in the formula.
    answer_c_expression = f"({s}-{one})({two}*{s}+{two}*{t}-{five})"
    answer_c_expression_simplified = "(s-1)(2s+2t-5)"


    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c_expression_simplified}")

solve_extremal_graph_problem()
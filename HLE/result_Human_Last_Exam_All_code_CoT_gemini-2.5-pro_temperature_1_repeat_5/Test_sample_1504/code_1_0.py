def solve_graph_theory_problem():
    """
    Solves the three-part graph theory problem and prints the formatted answer.
    """

    # Part (a): Determine if ex(n; G, K_{1,t}-ind) = Theta(n) when G is not a union of K2's.
    # The induced K_{1,t}-free property implies a maximum degree of at most t-1.
    # This immediately gives an upper bound of |E| <= n*(t-1)/2, so it's O(n).
    # For the lower bound, since G is not a union of K2's, it has a component with at least 3 vertices.
    # A graph H consisting of floor(n/2) disjoint edges (a matching) has Theta(n) edges.
    # This H has max degree 1, which is <= t-1 (since t>=2).
    # H is also G-free because all its components are K2.
    # Thus, the function is Theta(n).
    answer_a = "True"

    # Part (b): Determine if ex(n; sK2, K_{1,t}-ind) is independent of n.
    # An sK2-free graph has a maximum matching of size at most s-1.
    # An induced K_{1,t}-free graph has a maximum degree of at most t-1.
    # The vertices of any maximal matching form a vertex cover. Let M be a maximal matching.
    # The size of the vertex cover V(M) is 2*|M| <= 2*(s-1).
    # The number of edges is at most the size of this cover times the max degree.
    # |E| <= 2*(s-1) * (t-1).
    # Since the number of edges is bounded by a constant depending only on s and t,
    # the function is independent of n for large n.
    answer_b = "True"

    # Part (c): Provide an upper bound for ex(n; sK2, K_{1,t}-ind).
    # As derived for part (b), the number of edges is bounded by 2*(s-1)*(t-1).
    # This expression serves as a valid upper bound.
    answer_c = "2*(s-1)*(t-1)"

    # Print the final answer in the required format.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_graph_theory_problem()
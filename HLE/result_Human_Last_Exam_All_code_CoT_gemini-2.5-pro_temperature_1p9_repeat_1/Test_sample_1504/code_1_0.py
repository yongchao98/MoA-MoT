def solve_extremal_problem():
    """
    This function formulates and prints the solution to the graph theory problem.
    """

    # Part (a): True or false: If G is not a union of K_2's, then ex(n; G, K_{1,t}-ind) = Theta(n).
    # Analysis: This is False. Consider G = P_3 (path of 3 vertices). A graph is P_3-free if and only if
    # it is a disjoint union of cliques. Any clique is induced K_{1,t}-free for t >= 2.
    # To maximize edges in a P_3-free graph on n vertices, we can take a K_{n-1} and an isolated vertex.
    # This graph has O(n^2) edges, which is not Theta(n).
    answer_a = "False"

    # Part (b): If G ~ sK_2, is ex(n; sK_2, K_{1,t}-ind) independent of n?
    # Analysis: This is True. An sK_2-free graph H has matching number nu(H) <= s-1.
    # Using a maximum matching decomposition, V(H) = V(M) U I where M is the maximum matching
    # and I is an independent set. |V(M)| = 2*nu(H) <= 2*(s-1). The number of edges within V(M) is bounded.
    # The number of edges between V(M) and I is also bounded due to the induced K_{1,t}-free property.
    # Thus, the total number of edges is bounded by a constant depending on s and t, not n.
    answer_b = "True"

    # Part (c): For G ~ sK_2, express the upper bound for ex(n; sK_2, K_{1,t}-ind).
    # Analysis: From the logic in (b), the number of edges |E| is bounded by:
    # |E| <= C(2*k, 2) + 2*k*(t-1), where k = nu(H) <= s-1.
    # To maximize, we set k = s-1:
    # Bound = C(2*(s-1), 2) + 2*(s-1)*(t-1)
    #       = (s-1)*(2s-3) + 2*(s-1)*(t-1)
    #       = (s-1) * (2s-3 + 2t-2)
    #       = (s-1) * (2s+2t-5)
    # Note: In the final formula string, the asterisks are included for clarity.
    answer_c = "(s - 1) * (2*s + 2*t - 5)"

    final_answer_string = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"

    print(f"<<<{final_answer_string}>>>")

solve_extremal_problem()
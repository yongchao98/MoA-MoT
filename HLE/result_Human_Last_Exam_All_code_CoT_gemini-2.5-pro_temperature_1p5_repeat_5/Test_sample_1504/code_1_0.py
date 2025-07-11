def solve_extremal_graph_problem():
    """
    This function formulates and prints the answers to the graph theory problem.
    """
    # (a) True/False: If G is not a union of K2's, then ex(n; G, K1,t-ind) = Theta(n).
    answer_a = "True"

    # (b) True/False: If G ~ sK2 for s >= 2, ex(n; sK2, K1,t-ind) is independent of n.
    answer_b = "True"

    # (c) For G ~ sK2, express the upper bound for ex(n; sK2, K1,t-ind) in terms of s and t.
    # From our derivation, a valid upper bound is (s-1) * (2*s + t - 2).
    answer_c = "(s-1)*(2*s+t-2)"

    # Print the final combined answer.
    # The format required is (a) [True/False]; (b) [True/False]; (c) [expression].
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer)

solve_extremal_graph_problem()
<<< (a) True; (b) True; (c) (s-1)*(2*s+t-2) >>>
def solve_graph_theory_problem():
    """
    Solves the given Turan-type extremal graph theory problem
    and formats the answer as specified.
    """

    # Part (a): If G is not a union of K2's, then ex(n; G, K_{1,t}-ind) = Theta(n).
    # This is False. A counter-example with G=P3 shows the extremal function can be Theta(n^2).
    answer_a = "False"

    # Part (b): If G ~ sK2, ex(n; sK2, K_{1,t}-ind) is independent of n.
    # This is True. The number of edges can be bounded by a constant depending on s and t.
    answer_b = "True"

    # Part (c): Upper bound for ex(n; sK2, K_{1,t}-ind) in terms of s and t.
    # The derived upper bound is (s-1)(2s + 2t - 5).
    # The numbers in the final equation are 1, 2, 2, 5.
    num1 = 1
    num2_s = 2
    num2_t = 2
    num3 = 5
    expression_c = f"(s - {num1}) * ({num2_s}*s + {num2_t}*t - {num3})"

    # Format the final output string.
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {expression_c}"
    print(final_answer)

solve_graph_theory_problem()
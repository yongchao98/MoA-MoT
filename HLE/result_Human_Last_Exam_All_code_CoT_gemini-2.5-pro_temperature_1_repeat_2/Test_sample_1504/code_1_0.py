def solve_graph_theory_problem():
    """
    Solves the given graph theory problem and prints the formatted answer.
    """
    # Part (a): True or false: If G is not a union of K2's, then ex(n; G, K1,t-ind) = Theta(n).
    # Based on the reasoning, this is True.
    answer_a = "True"

    # Part (b): If G ~ sK2 for some s >= 2, is ex(n; sK2, K1,t-ind) independent of n?
    # Based on the reasoning, this is True.
    answer_b = "True"

    # Part (c): For G ~ sK2, express the upper bound for ex(n; sK2, K1,t-ind) in terms of s and t.
    # The derived upper bound is (s-1)(2s+2t-5).
    # The prompt requires outputting each number in the final equation.
    s_var = "s"
    t_var = "t"
    num_1 = 1
    num_2_s = 2
    num_2_t = 2
    num_5 = 5
    answer_c_expression = f"({s_var} - {num_1})({num_2_s}{s_var} + {num_2_t}{t_var} - {num_5})"

    # Format the final output string
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c_expression}"

    print(final_answer)

solve_graph_theory_problem()
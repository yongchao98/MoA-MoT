def solve_turan_problem():
    """
    This function formats the solution to the Turan-type extremal problem.
    The reasoning is provided in the text above.
    """
    
    # (a) True or false: If G is not a union of K2's, then ex(n; G, K1,t-ind) = Θ(n).
    answer_a = "True"
    
    # (b) If G ~ sK2 for some s ≥ 2, state whether the following is true or false:
    # ex(n; sK2, K1,t-ind) is independent of n.
    answer_b = "True"
    
    # (c) For G ~ sK2, express the upper bound for ex(n; sK2, K1,t-ind) in terms of s and t.
    # The expression is (s-1)(2s+2t-5).
    # We format it to show the numbers in the equation.
    s_var = "s"
    t_var = "t"
    num_1 = 1
    num_2_s = 2
    num_2_t = 2
    num_5 = 5
    answer_c = f"({s_var}-{num_1})*({num_2_s}*{s_var}+{num_2_t}*{t_var}-{num_5})"
    
    # The problem asks for the expression, not the expanded form with variables.
    # Let's provide the cleaner algebraic expression.
    answer_c_final = "(s-1)(2s+2t-5)"

    final_answer_string = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c_final}"
    
    print(final_answer_string)

solve_turan_problem()
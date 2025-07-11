def solve_extremal_problem():
    """
    This function formulates and prints the solution to the Turan-type extremal problem.
    """
    # Part (a): Based on the reasoning that the maximum degree is bounded,
    # and a linear construction exists.
    answer_a = "True"
    
    # Part (b): Based on the reasoning that the number of edges is bounded
    # by a constant dependent only on s and t.
    answer_b = "True"
    
    # Part (c): The expression for the upper bound.
    # We construct the formula string piece by piece as requested.
    s_minus_1_str = "(s - 1)"
    s_minus_2_str = "(s - 2)"
    t_minus_1_str = "(t - 1)"
    two_str = "2"
    
    # (s-1 choose 2) term, written as (s-1)*(s-2)/2
    term1_str = f"({s_minus_1_str} * {s_minus_2_str} / {two_str})"
    
    # (s-1)*(t-1) term
    term2_str = f"({s_minus_1_str} * {t_minus_1_str})"
    
    answer_c_expression = f"{term1_str} + {term2_str}"

    # Print the final result in the specified format
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c_expression}")

solve_extremal_problem()
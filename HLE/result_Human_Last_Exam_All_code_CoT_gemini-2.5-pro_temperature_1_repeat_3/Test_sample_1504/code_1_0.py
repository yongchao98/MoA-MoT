def solve_graph_theory_question():
    """
    This function provides the solution to the three-part graph theory question.
    """

    # Part (a): Based on the analysis, the statement is True.
    answer_a = "True"

    # Part (b): Based on the analysis, the statement is True.
    answer_b = "True"

    # Part (c): The expression for the upper bound.
    # The numbers in the expression are 1, 1, 1, 2, 2.
    s_minus_one = "s-1"
    t_minus_one = "t-1"
    s_minus_two = "s-2"
    two = "2"
    
    # We construct the expression string as requested.
    answer_c_expression = f"({s_minus_one})*({t_minus_one}) + ({s_minus_one})*({s_minus_two})/{two}"

    # Print the final answer in the specified format.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c_expression}")

solve_graph_theory_question()
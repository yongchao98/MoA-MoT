def solve_extremal_problem():
    """
    This function formulates and prints the answers to the extremal graph theory problem.
    """

    # Answer for part (a)
    answer_a = "False"

    # Answer for part (b)
    answer_b = "True"

    # Answer for part (c): The derived upper bound expression.
    # The format (s-1)(2s+2t-5) contains the numbers 1, 2, 2, and 5.
    expression_c = "(s-1) * (2*s + 2*t - 5)"

    # Print the final combined answer in the required format.
    final_answer_string = f"(a) {answer_a}; (b) {answer_b}; (c) {expression_c}"
    print(final_answer_string)
    
    # Do not remove the '<<<...>>>' marker. It is used for result extraction.
    print(f'<<<(a) False; (b) True; (c) (s-1)*(2*s + 2*t - 5)>>>')

solve_extremal_problem()
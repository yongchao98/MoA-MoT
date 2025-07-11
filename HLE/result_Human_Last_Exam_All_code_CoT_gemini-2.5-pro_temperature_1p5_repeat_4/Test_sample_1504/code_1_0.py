def solve_extremal_problem():
    """
    This function formulates and prints the solution to the Turan-type extremal problem.
    """

    # Part (a): True. The number of edges is of linear order Theta(n).
    answer_a = "True"

    # Part (b): False. The function is not independent of n for all n,
    # as its value changes for small n before becoming constant.
    answer_b = "False"

    # Part (c): The derived upper bound on the number of edges.
    # The expression is derived as Binomial(2*(s-1)*t, 2).
    # This expands to t*(s-1)*(2*t*(s-1) - 1).
    answer_c = "t*(s-1)*(2*t*(s-1) - 1)"

    # Format the final answer as requested
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"

    print(f"<<<{final_answer}>>>")

solve_extremal_problem()
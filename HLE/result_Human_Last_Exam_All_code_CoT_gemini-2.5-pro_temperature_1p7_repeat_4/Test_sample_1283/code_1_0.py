def solve_equation_roots():
    """
    This function calculates and prints the answers to the problem.
    """
    # Part (a): Determine the expression for the maximum number of solutions.
    # Based on the derivation, the maximum number of solutions is d_P + d_Q + 2.
    answer_a = "d_P + d_Q + 2"

    # Part (b): Calculate the maximum number for d_P = 3 and d_Q = 2.
    d_P = 3
    d_Q = 2

    # The formula from part (a) is applied.
    # The code below calculates the result based on the provided numbers.
    max_solutions = d_P + d_Q + 2

    # The required answer for (b) is the numerical result of the expression.
    answer_b = max_solutions

    print(f"(a) {answer_a}; (b) {answer_b}")

solve_equation_roots()
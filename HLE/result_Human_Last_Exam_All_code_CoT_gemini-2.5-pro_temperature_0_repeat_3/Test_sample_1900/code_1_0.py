def solve_complexity_question():
    """
    This function determines and prints the computational complexity for the described problem.

    The problem of finding a path through all n red balls is equivalent to the
    Hamiltonian Path problem. This problem is NP-complete, and the best-known
    algorithms have exponential time complexity.
    """

    # The base of the exponential term in the complexity O(base^n * n^exponent)
    base = 2

    # The exponent of the polynomial term in the complexity O(base^n * n^exponent)
    exponent = 2

    # Format the complexity string using the defined numbers
    complexity_string = f"O({base}^n * n^{exponent})"

    # The complexity is the same for both the decision problem (Question A)
    # and the search problem (Question B).
    answer_A = complexity_string
    answer_B = complexity_string

    # Print the final answer in the format "O(...); O(...)"
    print(f"{answer_A}; {answer_B}")

solve_complexity_question()
def solve_complexity_questions():
    """
    This function determines and prints the computational complexity for the given problems.
    The problem of finding a line through all balls is equivalent to the Hamiltonian Path problem
    in a graph. The properties of the red balls imply the graph is connected and locally connected.
    Hamiltonian Path is NP-complete even for this class of graphs, leading to exponential complexity.
    """

    # For Question A: Deciding if a path exists.
    # The complexity is exponential. A standard bound is O(n^2 * 2^n).
    poly_factor_A = 2
    exponential_base_A = 2
    complexity_A = f"O(n^{poly_factor_A} * {exponential_base_A}^n)"

    # For Question B: Finding the path if one exists.
    # The search problem is also NP-hard, with the same complexity as the decision problem.
    poly_factor_B = 2
    exponential_base_B = 2
    complexity_B = f"O(n^{poly_factor_B} * {exponential_base_B}^n)"

    # The final answer combines both complexities, separated by a semicolon.
    final_answer = f"{complexity_A}; {complexity_B}"
    print(final_answer)

solve_complexity_questions()
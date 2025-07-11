def solve_complexity_questions():
    """
    This function provides the computational complexity for the two questions based on the analysis above.
    The problem of finding a line through all n red balls is equivalent to the
    Hamiltonian Path Problem (HPP) in a graph with n vertices.

    HPP is an NP-complete problem, and even with the given constraints (grid graph structure,
    connectedness, and connected neighborhoods), it remains NP-complete.
    The best-known algorithms for solving HPP run in exponential time.

    A standard dynamic programming algorithm for HPP has a time complexity of O(2^n * n^2).
    """

    # Complexity for Question A: Deciding if a path exists.
    # The numbers in the equation are the base of the exponent and the power of the polynomial factor.
    base_A = 2
    poly_power_A = 2
    complexity_A = f"O({base_A}^n * n^{poly_power_A})"

    # Complexity for Question B: Finding the path.
    # The complexity is the same as the decision problem.
    base_B = 2
    poly_power_B = 2
    complexity_B = f"O({base_B}^n * n^{poly_power_B})"

    # The final answer is a semicolon-separated string of the two complexities.
    print(f"{complexity_A}; {complexity_B}")

solve_complexity_questions()
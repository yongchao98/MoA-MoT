def solve_complexity():
    """
    This function prints the computational complexity for the two questions.
    """
    # The complexity for Question A: Deciding if a path through all balls exists.
    # This is the Hamiltonian Path decision problem.
    # For a graph with n vertices, the complexity using dynamic programming is O(n^2 * 2^n).
    complexity_A = "O(n^2 * 2^n)"

    # The complexity for Question B: Finding such a path, given it exists.
    # This is the search version of the Hamiltonian Path problem.
    # The same dynamic programming algorithm can find the path with the same complexity.
    complexity_B = "O(n^2 * 2^n)"

    # The final answer should be in the format: complexity_A; complexity_B
    print(f"{complexity_A}; {complexity_B}")

solve_complexity()
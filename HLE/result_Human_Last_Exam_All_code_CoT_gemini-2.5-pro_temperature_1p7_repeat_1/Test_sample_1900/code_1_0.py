def solve_complexity_question():
    """
    This function determines and prints the computational complexity for the two questions.
    """
    # For Question A, the problem is to decide if a Hamiltonian path exists
    # in a graph that is guaranteed to be connected and locally connected.
    # A known theorem states that such a path always exists.
    # Therefore, the decision is trivial and takes constant time.
    # The number for the O(1) notation is 1.
    number_for_A = 1

    # For Question B, the problem is to find the Hamiltonian path that is
    # known to exist. For the special class of connected, locally connected graphs,
    # polynomial-time algorithms exist. A well-known complexity bound for
    # this search problem is O(n^3).
    # The number (power) in the final equation O(n^3) is 3.
    power_for_B = 3

    # Construct and print the final answer string.
    print(f"O({number_for_A}); O(n^{power_for_B})")

solve_complexity_question()
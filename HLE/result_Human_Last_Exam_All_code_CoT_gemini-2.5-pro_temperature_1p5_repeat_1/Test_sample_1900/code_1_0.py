def solve_complexity_game():
    """
    This function determines and prints the computational complexities for the described game.
    """
    # For Question A: Deciding if a line exists.
    # The given properties guarantee that the graph of red balls has a
    # Hamiltonian path. Therefore, the decision is always "yes" and can be
    # made without inspecting the graph. This is a constant time operation.
    complexity_A_value = 1
    answer_A = f"O({complexity_A_value})"

    # For Question B: Finding the line.
    # We need to find the Hamiltonian path. The properties of the graph,
    # combined with the fact that it's embedded on a grid (implying a
    # maximum degree of 8), allow for a polynomial-time algorithm.
    # A known cycle-expansion algorithm runs in O(n^2) time on such graphs.
    complexity_B_exponent = 2
    answer_B = f"O(n^{complexity_B_exponent})"

    # Combine the answers as requested, separated by a semicolon.
    final_answer = f"{answer_A}; {answer_B}"

    print(final_answer)

solve_complexity_game()
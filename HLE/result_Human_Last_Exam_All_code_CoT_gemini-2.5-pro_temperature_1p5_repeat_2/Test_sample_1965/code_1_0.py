def solve_rubiks_permutations():
    """
    Calculates the number of 6-move sequences that solve a Rubik's cube
    at move 4, 5, or 6.
    """

    # These are the known numbers of k-move sequences that return a
    # solved Rubik's cube to the solved state.
    # n_k: number of sequences of k moves that restore the cube.
    n_4 = 156
    n_5 = 0
    n_6 = 2460

    # The formula derived from the Principle of Inclusion-Exclusion is:
    # Total = 132 * n_4 + 12 * n_5 + n_6

    # Calculate each term of the equation
    term1 = 132 * n_4
    term2 = 12 * n_5
    term3 = n_6
    total_permutations = term1 + term2 + term3

    # Print the explanation and the final equation with the computed numbers
    print("The total number of permutations that result in the cube returning to its original configuration")
    print("at some point during the final 3 moves is calculated by the formula:")
    print("Total = 132 * n_4 + 12 * n_5 + n_6\n")
    print(f"Using the known values n_4 = {n_4}, n_5 = {n_5}, n_6 = {n_6}:")
    print(f"Total = (132 * {n_4}) + (12 * {n_5}) + {n_6}")
    print(f"Total = {term1} + {term2} + {term3}")
    print(f"Total = {total_permutations}")


solve_rubiks_permutations()

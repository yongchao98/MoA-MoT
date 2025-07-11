def solve_hat_puzzle():
    """
    Calculates the difference in guaranteed correct guesses between two scenarios.
    """
    # Total number of individuals
    num_individuals = 9

    # --- Calculation for N (Simultaneous Guessing) ---
    # The problem implies two independent binary properties for the hats.
    # For a single property, with 9 people, we can guarantee floor(9/2) = 4 correct guesses.
    # This is done by splitting the group into a team of 4 and a team of 5,
    # where each team bets on a different overall parity.
    #
    # For two properties, we need to find the guaranteed number of people correct for BOTH.
    # This involves finding the minimum overlap between the potentially correct groups.
    # Let's have two partitions into groups of 4 (A1, A2) and 5 (B1, B2).
    # The number of correct people can be |A1∩A2|, |A1∩B2|, |B1∩A2|, or |B1∩B2|.
    # We want to maximize the minimum of these values.
    # Let x = |A1 ∩ A2|. The intersection sizes are: x, 4-x, 4-x, and 1+x.
    # By testing x from 0 to 4, we find the maximum of min(x, 4-x, 1+x) occurs at x=2,
    # where the minimum value is 2.
    N = 2

    # --- Calculation for M (Sequential Guessing) ---
    # One person speaks first. They act as an informant.
    # They can see the other 8 hats. Their guess can encode 2 bits of information
    # (e.g., parity of property 1 for the 8 hats, parity of property 2 for the 8 hats).
    # The other 8 people hear this information, see the 7 other hats in their group,
    # and can deduce their own two hat properties with certainty.
    # This guarantees that n-1 people are correct.
    M = num_individuals - 1

    # --- Final Calculation ---
    difference = M - N

    print(f"The number of guaranteed correct guesses in the simultaneous scenario is N = {N}.")
    print(f"The number of guaranteed correct guesses in the sequential scenario is M = {M}.")
    print(f"The difference (M - N) is {M} - {N} = {difference}.")

solve_hat_puzzle()
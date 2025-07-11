def solve():
    """
    This function calculates and prints the answers to the user's questions.
    """
    # Part 1: Expected value of X_19
    n_19 = 19
    ex_19 = (n_19**2 - 1) / 2
    print(f"E[X_19] = {int(ex_19)}")

    # Part 2: Expected value of X_20
    # For even n, the expected time is infinite.
    ex_20 = "infinity"
    print(f"E[X_20] = {ex_20}")

    # Part 3: General formula for odd n
    # The formula is (n^2 - 1) / 2
    print("For odd n > 1, E[X_n] is given by the formula (n^2 - 1) / 2.")
    print("The numbers in the formula are: n, 2, 1, 2")

    # Part 4: General formula for expected number of times
    # The formula is 2 * (n - 13)
    print("For odd n > 30, the expected number of times the specified separation occurs is 2 * (n - 13).")
    print("The numbers in the formula are: 2, n, 13")
    
    # Part 5: Probability of game ending
    prob_one_ends = "Yes"
    print(f"Is it true that the game ends with probability one for every odd n > 1? {prob_one_ends}")

    # Final combined answer as requested by the user prompt
    final_answer = f"{int(ex_19)}, \u221e, (n^2-1)/2, 2(n-13), Yes"
    print("\nFormatted answer:")
    print(final_answer)


solve()
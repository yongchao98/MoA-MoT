import scipy.stats

def solve_coin_game():
    """
    Analyzes the coin game to determine if it's better to be the 1st or 2nd player.
    """
    # Game parameters
    num_1_euro = 136
    num_2_euro = 87
    total_coins = num_1_euro + num_2_euro

    # A line of 223 coins has 112 odd-indexed positions and 111 even-indexed positions.
    odd_positions = (total_coins + 1) // 2
    even_positions = total_coins // 2

    # The core of the strategy is that Player 1 (P1) can force the outcome where their score
    # is the sum of all coins at odd positions (Sum_odd), and Player 2 (P2) gets the sum
    # of all coins at even positions (Sum_even).

    # Let O2 be the number of 2-euro coins in the 112 odd positions.
    # The number of 1-euro coins in odd positions is (112 - O2).
    # P1's score is Sum_odd = 1 * (112 - O2) + 2 * O2 = 112 + O2.

    # Let E2 be the number of 2-euro coins in the 111 even positions.
    # The number of 1-euro coins in even positions is (111 - E2).
    # P2's score is Sum_even = 1 * (111 - E2) + 2 * E2 = 111 + E2.
    
    # We know that the total number of 2-euro coins is 87, so O2 + E2 = 87.

    print("--- Analysis of the Coin Game ---")
    print(f"There are {total_coins} coins in total: {num_1_euro} x 1-euro and {num_2_euro} x 2-euro.")
    print(f"P1 picks {odd_positions} coins and P2 picks {even_positions} coins.")
    print("\nWith optimal play, P1's score will be the sum of coins in odd positions, and P2's the sum of coins in even positions.")
    print("\nP1's score = (Value of coins in odd positions) = 1 * (112 - O2) + 2 * O2 = 112 + O2")
    print("P2's score = (Value of coins in even positions) = 1 * (111 - E2) + 2 * E2 = 111 + E2")
    print("where O2 and E2 are the number of 2-euro coins in odd and even positions, respectively.")

    print("\nP1 wins if P1's score > P2's score:")
    print("112 + O2 > 111 + E2")
    print("1 + O2 > E2")
    print("Substitute E2 = 87 - O2:")
    print("1 + O2 > 87 - O2")
    print("2 * O2 > 86")
    print("O2 > 43")

    print("\nThis means P1 wins if the number of 2-euro coins in the 112 odd positions is 44 or more.")
    print("We need to find the probability of this event for a random arrangement.")

    # The number of 2-euro coins in odd positions (O2) follows a hypergeometric distribution.
    # Population size (M): total_coins = 223
    # Number of successes in population (n): num_2_euro = 87
    # Sample size (N): odd_positions = 112
    M, n, N = total_coins, num_2_euro, odd_positions

    # The expected value of O2 gives a good indication of which outcome is more likely.
    expected_O2 = N * (n / M)
    print(f"\nThe expected number of 2-euro coins in odd positions is E[O2] = {N} * ({n}/{M}) = {expected_O2:.4f}.")
    print("Since the expected value (43.695) is greater than 43, it is more probable that O2 will be 44 or more.")

    # For a definitive answer, we calculate the probability P(O2 > 43).
    # This is 1 - P(O2 <= 43), which is 1 minus the cumulative distribution function (CDF) at k=43.
    prob_p1_wins = 1 - scipy.stats.hypergeom.cdf(k=43, M=M, n=n, N=N)

    print(f"\nThe precise probability of P1 winning is P(O2 > 43) = {prob_p1_wins:.4f}.")
    print("\nBecause this probability is greater than 0.5, you have a better chance of winning if you choose to be the 1st player.")

solve_coin_game()
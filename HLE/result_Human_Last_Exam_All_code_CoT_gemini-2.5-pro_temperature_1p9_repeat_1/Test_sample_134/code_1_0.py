import scipy.stats as stats

def solve_coin_game():
    """
    Analyzes the coin game to determine if it's better to be the 1st or 2nd player.
    """
    # 1. Game Parameters
    one_euro_coins = 136
    two_euro_coins = 87
    total_coins = one_euro_coins + two_euro_coins

    odd_positions = total_coins // 2 + 1
    even_positions = total_coins // 2

    # 2. Optimal Strategy Analysis & Winning Conditions
    # With optimal play, Player 2 (P2) wins if SumOdd < SumEven, where SumOdd is the
    # total value of coins in the odd_positions (1, 3, 5, ...) and SumEven is the total
    # value in the even_positions (2, 4, 6, ...).

    # Let n1_odd be the number of 1-euro coins in odd positions.
    # The number of 2-euro coins in odd positions is (odd_positions - n1_odd).
    # SumOdd = 1 * n1_odd + 2 * (odd_positions - n1_odd)

    # Let n1_even be the number of 1-euro coins in even positions.
    # n1_even = one_euro_coins - n1_odd
    # SumEven = 1 * n1_even + 2 * (even_positions - n1_even)

    # The inequality SumOdd < SumEven becomes:
    # 1*n1_odd + 2*(112 - n1_odd) < 1*(136 - n1_odd) + 2*(111 - (136 - n1_odd))
    # After simplification, this becomes:
    # 224 - n1_odd < 86 + n1_odd
    # 138 < 2 * n1_odd
    # n1_odd > 69

    # So, P2 has a guaranteed winning strategy if more than 69 one-euro coins
    # are randomly placed in the 112 odd positions. If n1_odd <= 69, P1 will
    # either win or draw (barring a few very unlikely edge cases).

    # 3. Probabilistic Calculation
    # We model n1_odd using the hypergeometric distribution.
    # Population size N: total positions
    # Number of successes in population M: number of odd positions
    # Sample size K: number of 1-euro coins being placed
    N = total_coins
    M = odd_positions
    K = one_euro_coins

    # The expected value for n1_odd
    expected_n1_odd = K * (M / N)

    # We want to find P(n1_odd > 69), which is 1 - P(n1_odd <= 69) = 1 - CDF(69)
    # CDF(k, N, M, K) gives P(X <= k)
    prob_p2_wins = 1 - stats.hypergeom.cdf(69, N, M, K)
    prob_p1_favored = 1 - prob_p2_wins

    # 4. Output the results
    print("Step 1: Game parameters")
    print(f"Total coins = {total_coins} ({one_euro_coins}x€1 + {two_euro_coins}x€2)")
    print(f"Number of odd positions = {odd_positions}")
    print(f"Number of even positions = {even_positions}")

    print("\nStep 2: Deriving Player 2's winning condition")
    print("Player 2 wins if SumOdd < SumEven.")
    print("This inequality can be expressed in terms of n1_odd (number of €1 coins on odd positions).")
    print("Equation for SumOdd: 1 * n1_odd + 2 * (112 - n1_odd) = 224 - n1_odd")
    print("Equation for SumEven: 1 * (136 - n1_odd) + 2 * (111 - (136 - n1_odd)) = 86 + n1_odd")
    print("Inequality: 224 - n1_odd < 86 + n1_odd")
    print("Solved for n1_odd: 138 < 2 * n1_odd  ==>  n1_odd > 69")
    print("Conclusion: Player 2 wins if 70 or more of the 136 €1 coins land on the 112 odd positions.")

    print("\nStep 3: Calculating probabilities")
    print("We use the hypergeometric distribution to find the probability of this event.")
    print(f"Parameters: N (total lots) = {N}, M (good lots) = {M}, K (draws) = {K}")
    print(f"The expected number of €1 coins on odd positions is E[n1_odd] = {K} * ({M}/{N}) = {expected_n1_odd:.4f}")
    print("\nThe probability that Player 2 wins (P(n1_odd > 69)) is:")
    print(f"P(P2 wins) = {prob_p2_wins:.4f}")
    print("\nThe probability that Player 1 wins or draws is:")
    print(f"P(P1 favored) = {prob_p1_favored:.4f}")
    
    print("\nStep 4: Conclusion")
    print("Since the probability of Player 1 having a winning position is much higher than for Player 2, you should choose to be the 1st player.")

solve_coin_game()
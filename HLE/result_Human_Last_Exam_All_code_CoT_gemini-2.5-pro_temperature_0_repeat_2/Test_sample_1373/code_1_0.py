def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    by making exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        float: The maximum probability of ending with 2*M GBP.
    """
    # The maximum possible amount of money we can have.
    # Starting with M, the max profit per trade is 12 (from Beta).
    max_money = M + n * 12

    # dp[i][j]: max probability of success with i trades remaining and j GBP.
    # Dimensions: (n + 1) for trades, (max_money + 1) for money.
    dp = [[0.0 for _ in range(max_money + 1)] for _ in range(n + 1)]

    # Base Case: i = 0 trades remaining.
    # Success is 1.0 if we have exactly 2*M, and 0.0 otherwise.
    # We must ensure 2*M is a reachable amount of money.
    if 0 <= 2 * M <= max_money:
        dp[0][2 * M] = 1.0

    # Fill the DP table bottom-up.
    # Iterate from i = 1 to n trades remaining.
    for i in range(1, n + 1):
        # Iterate through all possible money amounts j.
        for j in range(max_money + 1):
            # --- Calculate probability if we choose Strategy Alpha ---
            prob_alpha = 0.0
            if j >= 1:  # Check if we can afford the £1 fee.
                # Success: 60% chance, money becomes j+1.
                p_success_alpha = dp[i - 1][j + 1] if (j + 1) <= max_money else 0.0
                # Failure: 40% chance, money becomes j-1.
                p_failure_alpha = dp[i - 1][j - 1] if (j - 1) >= 0 else 0.0
                prob_alpha = 0.6 * p_success_alpha + 0.4 * p_failure_alpha

            # --- Calculate probability if we choose Strategy Beta ---
            prob_beta = 0.0
            if j >= 3:  # Check if we can afford the £3 fee.
                # Success: 20% chance, money becomes j+12.
                p_success_beta = dp[i - 1][j + 12] if (j + 12) <= max_money else 0.0
                # Failure: 80% chance, money becomes j-3.
                p_failure_beta = dp[i - 1][j - 3] if (j - 3) >= 0 else 0.0
                prob_beta = 0.2 * p_success_beta + 0.8 * p_failure_beta

            # The optimal strategy is to choose the one with the maximum probability.
            dp[i][j] = max(prob_alpha, prob_beta)

    # The final answer is the probability starting with M pounds and n trades.
    final_probability = dp[n][M]
    print(f"For an initial investment M={M} and n={n} trades:")
    print(f"The maximum probability of doubling the investment is: {final_probability}")

if __name__ == '__main__':
    # Example usage:
    # Let's use a case where a non-zero probability exists.
    # If M=13, n=2, we want to reach 26.
    # Path 1: Trade 1 (Alpha win, 60%) -> money=14. Trade 2 (Beta win, 20%) -> money=26. Prob = 0.6*0.2=0.12
    # Path 2: Trade 1 (Beta win, 20%) -> money=25. Trade 2 (Alpha win, 60%) -> money=26. Prob = 0.2*0.6=0.12
    # The optimal choice at the start gives a probability of 0.12.
    initial_investment_M = 13
    number_of_trades_n = 2
    solve_trading_probability(initial_investment_M, number_of_trades_n)
    
    print("-" * 20)
    
    # Another example
    # If M=12, n=1, we want to reach 24.
    # Only possible by winning with Beta (20% chance).
    initial_investment_M = 12
    number_of_trades_n = 1
    solve_trading_probability(initial_investment_M, number_of_trades_n)

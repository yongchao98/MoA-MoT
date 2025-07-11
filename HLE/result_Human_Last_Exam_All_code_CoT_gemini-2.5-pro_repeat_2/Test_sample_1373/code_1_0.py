def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        float: The maximum probability of reaching exactly 2*M GBP.
    """
    # The maximum possible capital is the initial capital plus the max profit from n trades.
    # Max profit per trade is from Strategy Beta: £15 return - £3 fee = £12.
    max_capital = M + 12 * n
    target_capital = 2 * M

    # dp[k][m] = max probability of success with k trades left and capital m.
    # Dimensions: (n+1) for trades, (max_capital+1) for money.
    dp = [[0.0 for _ in range(max_capital + 1)] for _ in range(n + 1)]

    # Base case: k = 0 trades remaining.
    # Success is 1.0 only if capital is exactly target_capital.
    if target_capital <= max_capital:
        dp[0][target_capital] = 1.0

    # Fill the DP table bottom-up.
    # Iterate through the number of trades remaining, from 1 to n.
    for k in range(1, n + 1):
        # Iterate through all possible capital amounts.
        for m in range(max_capital + 1):
            prob_alpha = 0.0
            prob_beta = 0.0

            # --- Calculate probability for Strategy Alpha ---
            if m >= 1:
                # Capital after success: m - 1 (fee) + 2 (return) = m + 1
                # Capital after failure: m - 1 (fee) + 0 (return) = m - 1
                p_succ = dp[k - 1][m + 1] if (m + 1) <= max_capital else 0.0
                p_fail = dp[k - 1][m - 1] if (m - 1) >= 0 else 0.0
                prob_alpha = 0.6 * p_succ + 0.4 * p_fail

            # --- Calculate probability for Strategy Beta ---
            if m >= 3:
                # Capital after success: m - 3 (fee) + 15 (return) = m + 12
                # Capital after failure: m - 3 (fee) + 0 (return) = m - 3
                p_succ = dp[k - 1][m + 12] if (m + 12) <= max_capital else 0.0
                p_fail = dp[k - 1][m - 3] if (m - 3) >= 0 else 0.0
                prob_beta = 0.2 * p_succ + 0.8 * p_fail
            
            # --- Determine optimal strategy for state (k, m) ---
            # The optimal probability is the max of the available options.
            dp[k][m] = max(prob_alpha, prob_beta)
            
    # The final answer is the probability at the start: n trades left, M capital.
    final_probability = dp[n][M]
    return final_probability

if __name__ == '__main__':
    # Example values for M and n
    # For M=13, n=2, the target is 26. A possible path to success is:
    # 1. Choose Alpha: 13 -> 14 (60% chance)
    # 2. Choose Beta: 14 -> 26 (20% chance)
    # Total probability for this sequence: 0.6 * 0.2 = 0.12
    # Another path:
    # 1. Choose Beta: 13 -> 25 (20% chance)
    # 2. Choose Alpha: 25 -> 26 (60% chance)
    # Total probability for this sequence: 0.2 * 0.6 = 0.12
    # The optimal strategy will pick one of these paths, resulting in a 0.12 probability.
    initial_investment = 13
    num_trades = 2

    result = solve_trading_probability(initial_investment, num_trades)
    
    print(f"Initial Investment (M): {initial_investment}")
    print(f"Number of Trades (n): {num_trades}")
    print(f"Target Capital (2M): {2 * initial_investment}")
    print(f"The maximum probability of success is: {result}")
    
    print("\n--- Another Example ---")
    initial_investment_2 = 9
    num_trades_2 = 2
    result_2 = solve_trading_probability(initial_investment_2, num_trades_2)
    print(f"Initial Investment (M): {initial_investment_2}")
    print(f"Number of Trades (n): {num_trades_2}")
    print(f"Target Capital (2M): {2 * initial_investment_2}")
    print(f"The maximum probability of success is: {result_2}")

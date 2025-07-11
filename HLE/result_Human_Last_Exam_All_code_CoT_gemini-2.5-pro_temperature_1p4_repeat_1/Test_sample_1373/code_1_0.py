def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        float: The maximum probability of success.
    """
    target_capital = 2 * M

    # The maximum profit in a single trade is +12 GBP (from Beta: 15 return - 3 fee).
    # A safe upper bound for capital is the initial amount plus max profit from all trades.
    max_capital = M + n * 12

    # We initialize the DP table.
    # dp[i][j] = max probability of success, given we have completed `i` trades and have `j` capital.
    # Dimensions: (n + 1) trades x (max_capital + 1) capital values.
    dp = [[0.0 for _ in range(max_capital + 1)] for _ in range(n + 1)]

    # Base case: At the end of n trades (i=n).
    # Success is only achieved if capital is exactly the target.
    # We must check if the target_capital is within our calculated bounds.
    if target_capital <= max_capital:
        dp[n][target_capital] = 1.0

    # Iterate backwards from the (n-1)-th trade down to the 0-th trade.
    for i in range(n - 1, -1, -1):
        # Iterate over all possible capital amounts at this stage.
        for j in range(max_capital + 1):
            
            # --- Possibility 1: Strategy Alpha ---
            prob_alpha = 0.0
            if j >= 1:  # Must have at least £1 for the fee.
                # Probability of future success if Alpha trade is successful (capital -> j+1)
                capital_s_a = j + 1
                prob_from_success_a = dp[i + 1][capital_s_a] if capital_s_a <= max_capital else 0.0

                # Probability of future success if Alpha trade fails (capital -> j-1)
                capital_f_a = j - 1
                prob_from_fail_a = dp[i + 1][capital_f_a] if capital_f_a >= 0 else 0.0
                
                prob_alpha = 0.6 * prob_from_success_a + 0.4 * prob_from_fail_a

            # --- Possibility 2: Strategy Beta ---
            prob_beta = 0.0
            if j >= 3:  # Must have at least £3 for the fee.
                # Probability of future success if Beta trade is successful (capital -> j+12)
                capital_s_b = j + 12
                prob_from_success_b = dp[i + 1][capital_s_b] if capital_s_b <= max_capital else 0.0

                # Probability of future success if Beta trade fails (capital -> j-3)
                capital_f_b = j - 3
                prob_from_fail_b = dp[i + 1][capital_f_b] if capital_f_b >= 0 else 0.0

                prob_beta = 0.2 * prob_from_success_b + 0.8 * prob_from_fail_b

            # The optimal strategy for state (i, j) is the one yielding the higher probability.
            dp[i][j] = max(prob_alpha, prob_beta)

    # The final answer is the probability at the start (0 trades done, M capital).
    return dp[0][M]

# --- Main Execution Block ---
# You can change these values to test different scenarios.
# Note: large values for M and n will increase computation time and memory usage.
initial_investment_M = 25
num_trades_n = 10

# Calculate the probability
optimal_probability = solve_trading_problem(initial_investment_M, num_trades_n)

# Print the final numbers and the resulting probability
print(f"Initial Investment (M): {initial_investment_M} GBP")
print(f"Number of Trades (n): {num_trades_n}")
print(f"Target Investment (2M): {2 * initial_investment_M} GBP")
print(f"The maximum probability of reaching the target is: {optimal_probability}")
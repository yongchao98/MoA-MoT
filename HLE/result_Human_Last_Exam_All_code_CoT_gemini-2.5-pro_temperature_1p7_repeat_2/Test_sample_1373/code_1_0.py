def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    by executing exactly n trades.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.
    """
    if M < 0 or n < 0:
        print("Initial investment and number of trades must be non-negative.")
        return

    target_money = 2 * M

    # The maximum possible net gain from one trade is +12 (from Beta strategy).
    # Starting with M, after n trades, the maximum possible money is M + 12*n.
    # This determines the maximum index needed for our money state.
    max_m = M + 12 * n

    # T[l][m] = max probability of success, given l trades are done and we have m GBP.
    # Dimensions: (n+1) for trades, (max_m+1) for money.
    T = [[0.0] * (max_m + 1) for _ in range(n + 1)]

    # Base case: After n trades are complete (l = n).
    # Success is defined as having exactly target_money.
    if target_money <= max_m:
        T[n][target_money] = 1.0

    # Fill the DP table by iterating backwards from l = n-1 down to 0.
    for l in range(n - 1, -1, -1):
        # Iterate through all possible money amounts at this stage.
        # The max money possible after l trades is M + 12*l, but we iterate over
        # the full range for simplicity.
        for m in range(max_m + 1):
            # Calculate probability of success if we choose Strategy Alpha
            # This is only possible if we can afford the £1 fee.
            prob_alpha = 0.0
            if m >= 1:
                # Get probabilities from the state after this trade (l+1)
                # Handle boundary conditions for money (cannot exceed max_m)
                p_success = T[l+1][m + 1] if (m + 1) <= max_m else 0.0
                p_failure = T[l+1][m - 1] if (m - 1) >= 0 else 0.0
                prob_alpha = 0.6 * p_success + 0.4 * p_failure

            # Calculate probability of success if we choose Strategy Beta
            # This is only possible if we can afford the £3 fee.
            prob_beta = 0.0
            if m >= 3:
                # Get probabilities from the state after this trade (l+1)
                p_success = T[l+1][m + 12] if (m + 12) <= max_m else 0.0
                p_failure = T[l+1][m - 3] if (m - 3) >= 0 else 0.0
                prob_beta = 0.2 * p_success + 0.8 * p_failure

            # Determine T[l][m] based on the optimal choice for the next trade.
            if m < 1:
                T[l][m] = 0.0  # Cannot afford any trade, so we can't complete n trades.
            elif m < 3:
                T[l][m] = prob_alpha  # Can only afford Alpha.
            else:
                T[l][m] = max(prob_alpha, prob_beta)  # Choose the better strategy.

    # The final answer is the probability at the start: 0 trades done, M capital.
    final_probability = T[0][M]

    print(f"Initial Investment (M): {M} GBP")
    print(f"Number of Trades (n): {n}")
    print(f"Target Investment: {target_money} GBP")
    print("-" * 40)
    print(f"The maximum probability of reaching exactly {target_money} GBP is: {final_probability:.6f}")


if __name__ == '__main__':
    # --- Example Usage ---
    # You can change these values to test different scenarios.
    initial_investment_M = 25
    number_of_trades_n = 10
    solve_trading_problem(initial_investment_M, number_of_trades_n)

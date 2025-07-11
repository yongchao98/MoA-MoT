def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    by executing exactly n trades.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.
    """
    if M <= 0 or n <= 0:
        print("Initial investment and number of trades must be positive.")
        return

    # The maximum possible profit in a single trade is £12 (from Beta).
    # So, the maximum possible money is M + n * 12.
    max_money = M + n * 12
    
    # DP table: dp[l][m] corresponds to T[l, m] from the problem description.
    # Dimensions: (n + 1) trades by (max_money + 1) pounds.
    dp = [[0.0 for _ in range(max_money + 1)] for _ in range(n + 1)]

    # Base Case: After n trades are complete.
    target_money = 2 * M
    if target_money <= max_money:
        # We succeed only if we have exactly 2M.
        dp[n][target_money] = 1.0

    # DP Calculation: Iterate backwards from the second-to-last trade to the first.
    for l in range(n - 1, -1, -1):
        # Iterate over all possible money amounts at step l.
        # Max money possible after l trades is M + l*12. We can loop to max_money for simplicity.
        for m in range(max_money + 1):
            
            # --- Probability of success if we choose Strategy Alpha ---
            p_alpha = 0.0
            if m >= 1:  # Must have at least £1 to invest.
                prob_succ_alpha = 0.6 * dp[l + 1][m + 1]
                prob_fail_alpha = 0.4 * (dp[l + 1][m - 1] if m - 1 >= 0 else 0)
                p_alpha = prob_succ_alpha + prob_fail_alpha

            # --- Probability of success if we choose Strategy Beta ---
            p_beta = 0.0
            if m >= 3:  # Must have at least £3 to invest.
                prob_succ_beta = 0.2 * dp[l + 1][m + 12]
                prob_fail_beta = 0.8 * (dp[l + 1][m - 3] if m - 3 >= 0 else 0)
                p_beta = prob_succ_beta + prob_fail_beta
            
            # The optimal strategy maximizes the probability of success.
            dp[l][m] = max(p_alpha, p_beta)
            
    # The final answer is the probability at the start (0 trades, M money).
    final_prob = dp[0][M]

    # --- Output the final equation as requested ---
    print(f"Let T[l, m] be the max probability of success with £m after l trades.")
    print(f"To find the answer, we compute T[0, {M}].\n")
    print(f"The decision at the start (l=0, m={M}) involves comparing two strategies:")

    # Calculate and display the Alpha option
    if M >= 1:
        p_alpha_final = 0.6 * dp[1][M + 1] + 0.4 * (dp[1][M - 1] if M > 0 else 0)
        print(f"  Choice 1: Strategy Alpha")
        print(f"    P(Success) = 0.60 * T[1, {M + 1}] + 0.40 * T[1, {M - 1}]")
        print(f"               = 0.60 * {dp[1][M+1]:.4f} + 0.40 * {(dp[1][M-1] if M > 0 else 0):.4f}")
        print(f"               = {p_alpha_final:.6f}")
    else:
        print("  Choice 1: Strategy Alpha is not affordable.")
        p_alpha_final = 0.0

    # Calculate and display the Beta option
    if M >= 3:
        p_beta_final = 0.2 * dp[1][M + 12] + 0.8 * (dp[1][M - 3] if M >= 3 else 0)
        print(f"  Choice 2: Strategy Beta")
        print(f"    P(Success) = 0.20 * T[1, {M + 12}] + 0.80 * T[1, {M - 3}]")
        print(f"               = 0.20 * {dp[1][M+12]:.4f} + 0.80 * {dp[1][M-3]:.4f}")
        print(f"               = {p_beta_final:.6f}")
    else:
        print("  Choice 2: Strategy Beta is not affordable.")
        p_beta_final = 0.0

    print(f"\nT[0, {M}] = max(P_alpha, P_beta) = max({p_alpha_final:.6f}, {p_beta_final:.6f})")
    print("\n---------------------------------------------------------")
    print(f"The maximum probability of success is: {final_prob:.6f}")
    print("---------------------------------------------------------")


if __name__ == '__main__':
    # You can change these values to test different scenarios
    initial_investment_M = 25
    number_of_trades_n = 10
    
    print(f"Calculating for Initial Investment M = £{initial_investment_M} and Number of Trades n = {number_of_trades_n}")
    print("-" * 30)
    solve_trading_probability(initial_investment_M, number_of_trades_n)
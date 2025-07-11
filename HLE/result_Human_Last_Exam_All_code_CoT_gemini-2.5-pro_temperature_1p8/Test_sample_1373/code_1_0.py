def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    Args:
        M (int): Initial investment in GBP.
        n (int): The exact number of trades to be executed.
    """
    # Define the target amount of money
    target_money = 2 * M

    # The maximum profit from a single trade is from a successful Beta strategy: £15 - £3 = £12.
    # So, the maximum possible capital after n trades is M + n * 12.
    max_capital = M + n * 12

    # Initialize DP table: dp[i][j] = max probability with i trades left and j capital.
    # Dimensions: (n+1) rows for trades [0...n], (max_capital+1) cols for money [0...max_capital]
    dp = [[0.0 for _ in range(max_capital + 1)] for _ in range(n + 1)]

    # Base case: When 0 trades are left (i=0).
    # Success is only possible if the capital is exactly the target amount.
    if target_money <= max_capital:
        dp[0][target_money] = 1.0

    # Fill the DP table bottom-up, from i=1 to n trades remaining.
    for i in range(1, n + 1):  # i = number of trades remaining
        for j in range(max_capital + 1):  # j = current capital

            # --- Option 1: Strategy Alpha ---
            prob_alpha = 0.0
            if j >= 1:  # Must have at least £1 to afford the fee
                # Probability of success if the trade returns £2 (net profit £1)
                p_success_from_alpha_win = 0
                if j + 1 <= max_capital:
                    p_success_from_alpha_win = dp[i-1][j+1]
                
                # Probability of success if the trade returns £0 (net loss £1)
                p_success_from_alpha_loss = 0
                if j - 1 >= 0:
                    p_success_from_alpha_loss = dp[i-1][j-1]
                
                prob_alpha = 0.60 * p_success_from_alpha_win + 0.40 * p_success_from_alpha_loss

            # --- Option 2: Strategy Beta ---
            prob_beta = 0.0
            if j >= 3:  # Must have at least £3 to afford the fee
                # Probability of success if the trade returns £15 (net profit £12)
                p_success_from_beta_win = 0
                if j + 12 <= max_capital:
                    p_success_from_beta_win = dp[i-1][j+12]

                # Probability of success if the trade returns £0 (net loss £3)
                p_success_from_beta_loss = 0
                if j - 3 >= 0:
                    p_success_from_beta_loss = dp[i-1][j-3]

                prob_beta = 0.20 * p_success_from_beta_win + 0.80 * p_success_from_beta_loss

            # The optimal choice for dp[i][j] is the max of the two strategies
            dp[i][j] = max(prob_alpha, prob_beta)

    # The final answer is the probability starting with M capital and n trades.
    final_probability = dp[n][M]

    # Print the problem setup and the result
    print(f"Problem Parameters:")
    print(f"Initial Investment (M): {M} GBP")
    print(f"Number of Trades (n): {n}")
    print(f"Target Investment: {target_money} GBP")
    
    print("\nDP Recurrence Relation:")
    print("T[m, l] = max(prob_alpha, prob_beta), where l is trades completed.")
    print("  prob_alpha = 0.60 * T[m+1, l+1] + 0.40 * T[m-1, l+1]")
    print("  prob_beta  = 0.20 * T[m+12, l+1] + 0.80 * T[m-3, l+1]")

    print(f"\nResult:")
    print(f"The maximum probability of reaching exactly {target_money} GBP is: {final_probability:.8f}")


if __name__ == '__main__':
    # --- Example Case ---
    # You can change these values to test different scenarios.
    initial_investment_M = 25
    number_of_trades_n = 10

    solve_trading_probability(initial_investment_M, number_of_trades_n)
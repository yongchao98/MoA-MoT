def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    by making exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be made.
    """
    
    # --- State Definition ---
    # Let dp[i][j] be the maximum probability of success (reaching 2*M at the end)
    # given that we have 'i' trades remaining and a current capital of 'j' GBP.

    # --- DP Table Size ---
    # A safe upper bound for money to track is determined by the maximum capital
    # from which it's still possible to reach the target 2*M. A loose but safe
    # bound is 2*M + 3*n. We add a buffer of 15 for index safety, as we need to
    # access indices like j+12 during calculations.
    MAX_MONEY = 2 * M + 3 * n + 15
    TARGET_MONEY = 2 * M

    # Initialize DP table with probabilities of 0.
    # Dimensions: (n+1) trades x MAX_MONEY capital states.
    dp = [[0.0] * MAX_MONEY for _ in range(n + 1)]

    # --- Base Case ---
    # With 0 trades remaining (i=0), the probability of success is 1.0 if we have
    # exactly the target capital, and 0.0 otherwise.
    if TARGET_MONEY >= 0 and TARGET_MONEY < MAX_MONEY:
        dp[0][TARGET_MONEY] = 1.0

    # --- Recurrence Relation ---
    # Fill the table bottom-up, from i=1 to n trades remaining.
    for i in range(1, n + 1):
        # Iterate over all relevant capital states 'j'. The loop range must ensure
        # that j+12 remains a valid index.
        for j in range(MAX_MONEY - 13):
            # A trade is only possible if capital is positive.
            if j <= 0:
                continue

            # Option 1: Strategy Alpha (fee: 1 GBP), available if j >= 1
            prob_from_alpha = 0.0
            # Net profit is +1 (success) or -1 (failure).
            prob_from_alpha = (0.60 * dp[i - 1][j + 1]) + (0.40 * dp[i - 1][j - 1])

            # Option 2: Strategy Beta (fee: 3 GBP), available if j >= 3
            prob_from_beta = 0.0
            if j >= 3:
                # Net profit is +12 (success) or -3 (failure).
                prob_from_beta = (0.20 * dp[i - 1][j + 12]) + (0.80 * dp[i - 1][j - 3])
            
            # The optimal strategy at state (i, j) is to choose the action
            # that maximizes the probability of success. If a strategy is not
            # available, its probability is 0.
            dp[i][j] = max(prob_from_alpha, prob_from_beta)

    # The final answer is the probability at the initial state: n trades, M capital.
    final_probability = dp[n][M]

    # --- Output ---
    print("This script solves a dynamic programming problem for a trading scenario.")
    print("\n--- Problem Setup ---")
    print(f"Initial Investment (M): {M} GBP")
    print(f"Number of Trades (n): {n}")
    print(f"Target Investment (2M): {2 * M} GBP")

    print("\n--- Final Equation (Recurrence Relation) ---")
    print("Let T[m, l] be the max probability of success with 'm' GBP and 'l' trades remaining.")
    print("The optimal strategy follows:")
    print("T[m, l] = max(Prob_from_Alpha, Prob_from_Beta)")
    print("\nWhere the numbers in the final equation for each strategy are:")
    print("Prob_from_Alpha = (0.60 * T[m + 1, l - 1]) + (0.40 * T[m - 1, l - 1])")
    print("Prob_from_Beta  = (0.20 * T[m + 12, l - 1]) + (0.80 * T[m - 3, l - 1])")

    print("\n--- Result ---")
    print(f"The maximum probability of achieving the target is: {final_probability:.6f}")


# --- Main execution ---
if __name__ == '__main__':
    # Example values for M and n. You can change these to test other scenarios.
    INITIAL_INVESTMENT_M = 25
    NUMBER_OF_TRADES_N = 10
    
    solve_trading_problem(INITIAL_INVESTMENT_M, NUMBER_OF_TRADES_N)
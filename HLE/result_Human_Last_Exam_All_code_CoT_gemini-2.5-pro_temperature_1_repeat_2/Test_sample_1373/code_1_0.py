import sys

def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.
    """
    # Increase recursion limit for potentially deep state dependencies, though our iterative approach avoids this.
    # sys.setrecursionlimit(max(sys.getrecursionlimit(), n + 5))

    # The target capital we want to reach.
    target_capital = 2 * M
    
    # The maximum possible capital is M + n * 12 (if all n trades are successful Beta trades).
    # This determines the size of our DP table's capital dimension.
    max_capital = M + n * 12
    
    # dp[i][j] = max probability of success with i trades left and j capital.
    # Initialize the DP table with 0.0.
    dp = [[0.0 for _ in range(max_capital + 1)] for _ in range(n + 1)]
    
    # policy[i][j] stores the optimal move ('A' for Alpha, 'B' for Beta, 'N' for None).
    policy = [['N' for _ in range(max_capital + 1)] for _ in range(n + 1)]

    # --- Base Case: i = 0 trades left ---
    # The probability of success is 1.0 if our capital is exactly the target,
    # and 0.0 otherwise.
    if target_capital <= max_capital:
        dp[0][target_capital] = 1.0

    # --- Iterative Calculation ---
    # We fill the DP table starting from i=1 up to n trades remaining.
    for i in range(1, n + 1):
        # Iterate through all possible capital amounts j.
        for j in range(max_capital + 1):
            
            # --- Evaluate Strategy Alpha ---
            # This strategy has a £1 fee.
            prob_alpha = -1.0  # Sentinel value for an impossible move.
            if j >= 1:
                # 60% chance of returning £2 (net profit £1) -> capital becomes j+1
                p_success_alpha = dp[i-1][j+1] if (j + 1) <= max_capital else 0.0
                # 40% chance of returning £0 (net loss £1) -> capital becomes j-1
                p_fail_alpha = dp[i-1][j-1] if (j - 1) >= 0 else 0.0
                prob_alpha = 0.60 * p_success_alpha + 0.40 * p_fail_alpha

            # --- Evaluate Strategy Beta ---
            # This strategy has a £3 fee.
            prob_beta = -1.0 # Sentinel value for an impossible move.
            if j >= 3:
                # 20% chance of returning £15 (net profit £12) -> capital becomes j+12
                p_success_beta = dp[i-1][j+12] if (j + 12) <= max_capital else 0.0
                # 80% chance of returning £0 (net loss £3) -> capital becomes j-3
                p_fail_beta = dp[i-1][j-3] if (j - 3) >= 0 else 0.0
                prob_beta = 0.20 * p_success_beta + 0.80 * p_fail_beta

            # --- Determine Optimal Strategy ---
            # Choose the move that maximizes the probability of success.
            if prob_alpha > prob_beta:
                dp[i][j] = prob_alpha
                policy[i][j] = 'A'
            elif prob_beta >= 0: # Check if Beta was a valid option
                dp[i][j] = prob_beta
                policy[i][j] = 'B'
            # If neither move is possible, dp[i][j] remains 0.0 and policy remains 'N'.

    # The final answer is the probability at the starting state (n trades, M capital).
    final_probability = dp[n][M]
    
    print(f"Initial Capital M = £{M}, Number of Trades n = {n}")
    print(f"Target Capital = £{target_capital}")
    print("-" * 40)

    print(f"Optimal probability of success: {final_probability:.6f}")

    if final_probability == 0.0:
        print("It is not possible to guarantee reaching the target capital.")
    else:
        # Show the calculation for the first move based on the policy table.
        first_move = policy[n][M]
        print(f"The optimal first move is Strategy {first_move} ('A' for Alpha, 'B' for Beta).")
        
        # Print the equation for the first step's probability calculation.
        if first_move == 'A':
            # Probability after one successful Alpha trade
            p1 = dp[n-1][M+1] if (M + 1) <= max_capital else 0.0
            # Probability after one failed Alpha trade
            p2 = dp[n-1][M-1] if (M - 1) >= 0 else 0.0
            # The calculation itself
            calc = 0.60 * p1 + 0.40 * p2
            print("Final Probability = 0.60 * P(n-1, M+1) + 0.40 * P(n-1, M-1)")
            print(f"Final Probability = 0.60 * {p1:.6f} + 0.40 * {p2:.6f} = {calc:.6f}")
        elif first_move == 'B':
            # Probability after one successful Beta trade
            p1 = dp[n-1][M+12] if (M + 12) <= max_capital else 0.0
            # Probability after one failed Beta trade
            p2 = dp[n-1][M-3] if (M - 3) >= 0 else 0.0
            # The calculation itself
            calc = 0.20 * p1 + 0.80 * p2
            print("Final Probability = 0.20 * P(n-1, M+12) + 0.80 * P(n-1, M-3)")
            print(f"Final Probability = 0.20 * {p1:.6f} + 0.80 * {p2:.6f} = {calc:.6f}")


if __name__ == '__main__':
    # --- Set initial parameters here ---
    # M: Initial capital in GBP
    initial_capital_M = 25
    # n: Exact number of trades to execute
    num_trades_n = 5
    
    solve_trading_problem(initial_capital_M, num_trades_n)

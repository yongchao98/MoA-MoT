def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        None. Prints the result and the calculation for the final step.
    """
    if M <= 0 or n < 0:
        print("Initial investment and number of trades must be positive.")
        return

    # Step 1: Define the DP state space.
    # money_range must be large enough for the target 2*M and the max possible holding M + 12*n.
    # We add a buffer of 13 to safely access indices like m+12 without extra boundary checks.
    money_range = max(2 * M, M + 12 * n) + 13

    # dp[l][m]: max probability of success with l trades left and m money.
    dp = [[0.0 for _ in range(money_range)] for _ in range(n + 1)]

    # Step 2: Set the base case (l=0 trades left).
    # Success is 1.0 if we have exactly 2*M, otherwise 0.
    if 2 * M < money_range:
        dp[0][2 * M] = 1.0

    # Step 3: Fill the DP table iteratively.
    for l in range(1, n + 1):
        # Determine a practical upper bound for money 'm' at this stage 'l'
        # to avoid iterating over unnecessarily large empty parts of the DP table.
        # Max money you could have started with to reach any state at l-1 is M + 12*(n-(l-1)).
        max_m_to_check = M + 12 * (n - (l - 1)) + 5 

        for m in range(1, min(money_range - 12, max_m_to_check)):
            
            prob_alpha = 0.0
            prob_beta = 0.0

            # Strategy Alpha: Fee £1. Net gain +£1 (60%), -£1 (40%). Requires m >= 1.
            val_alpha_success = dp[l-1][m + 1]
            val_alpha_fail = dp[l-1][m - 1]
            prob_alpha = 0.6 * val_alpha_success + 0.4 * val_alpha_fail

            # Strategy Beta: Fee £3. Net gain +£12 (20%), -£3 (80%). Requires m >= 3.
            if m >= 3:
                val_beta_success = dp[l-1][m + 12]
                val_beta_fail = dp[l-1][m - 3]
                prob_beta = 0.2 * val_beta_success + 0.8 * val_beta_fail
            
            # Optimal choice is the max of the probabilities of affordable strategies.
            dp[l][m] = max(prob_alpha, prob_beta)

    # Step 4: The result is the state for the start of the problem.
    final_result = dp[n][M]
    print(f"Initial Investment (M): {M} GBP")
    print(f"Number of Trades (n): {n}")
    print(f"Target Investment (2M): {2*M} GBP")
    print(f"Maximum probability of reaching exactly {2*M} GBP is: {final_result}")
    
    # As requested, output the numbers in the "final equation", which corresponds
    # to the calculation for dp[n][M] (the first trade).
    if n > 0:
        final_p_alpha = 0.0
        final_p_beta = 0.0

        print("\n--- Final Step Calculation ---")
        print("To find the answer, the algorithm computes the optimal choice for the first trade, starting with")
        print(f"M={M} and n={n} trades remaining. This relies on the pre-computed probabilities for n-1 trades.")
        
        if M >= 1:
            dp_val_alpha_succ = dp[n - 1][M + 1]
            dp_val_alpha_fail = dp[n - 1][M - 1]
            final_p_alpha = 0.6 * dp_val_alpha_succ + 0.4 * dp_val_alpha_fail
            print(f"\nIf choosing Alpha:")
            print(f"  Success changes money to £{M+1}, fail to £{M-1}. Remaining trades: {n-1}.")
            print(f"  Prob(Success from £{M+1}, {n-1} trades) = {dp_val_alpha_succ}")
            print(f"  Prob(Success from £{M-1}, {n-1} trades) = {dp_val_alpha_fail}")
            print(f"  Resulting Total Probability = 0.6 * {dp_val_alpha_succ} + 0.4 * {dp_val_alpha_fail} = {final_p_alpha}")

        if M >= 3:
            dp_val_beta_succ = dp[n-1][M+12]
            dp_val_beta_fail = dp[n-1][M-3]
            final_p_beta = 0.2 * dp_val_beta_succ + 0.8 * dp_val_beta_fail
            print(f"\nIf choosing Beta:")
            print(f"  Success changes money to £{M+12}, fail to £{M-3}. Remaining trades: {n-1}.")
            print(f"  Prob(Success from £{M+12}, {n-1} trades) = {dp_val_beta_succ}")
            print(f"  Prob(Success from £{M-3}, {n-1} trades) = {dp_val_beta_fail}")
            print(f"  Resulting Total Probability = 0.2 * {dp_val_beta_succ} + 0.8 * {dp_val_beta_fail} = {final_p_beta}")
        
        print("\nThe algorithm chooses the strategy with the maximum resulting probability to determine the value.")

# Example Usage: You can change these values to test different scenarios.
# Note: For very large M or n, this will be slow and memory-intensive.
initial_investment_M = 25
num_trades_n = 10
solve_trading_problem(initial_investment_M, num_trades_n)
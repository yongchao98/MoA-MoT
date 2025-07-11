import sys

def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        float: The maximum probability of reaching exactly 2*M GBP.
    """
    if M <= 0 or n < 0:
        return 0.0

    # The maximum possible capital we can have is if we start with M
    # and every trade is a successful Beta trade (+£12).
    # We add a small buffer for safety.
    MAX_CAPITAL = M + 12 * n + 5
    TARGET_CAPITAL = 2 * M

    if TARGET_CAPITAL > MAX_CAPITAL:
        # It's impossible to reach the target if it's beyond the max possible capital
        # This check is mostly for sanity, the DP logic would handle it anyway.
        print(f"Target capital {TARGET_CAPITAL} is unreachable.")
        return 0.0

    # dp[i][j]: max probability of success with i trades remaining and j pounds.
    # Dimensions: (n+1) rows for trades, (MAX_CAPITAL+1) columns for money.
    dp = [[0.0] * (MAX_CAPITAL + 1) for _ in range(n + 1)]

    # Base Case: i = 0 trades remaining.
    # Success is 1.0 only if capital is exactly the target, 0.0 otherwise.
    if TARGET_CAPITAL <= MAX_CAPITAL:
        dp[0][TARGET_CAPITAL] = 1.0

    # Fill the DP table bottom-up.
    # Iterate through the number of trades remaining, from 1 to n.
    for i in range(1, n + 1):
        # Iterate through each possible capital amount j.
        for j in range(MAX_CAPITAL + 1):
            
            # --- Calculate probability of success if we choose Strategy Alpha ---
            prob_alpha = 0.0
            # We must be able to afford the £1 fee.
            if j >= 1:
                # The final equation for the probability of success from this state (i, j)
                # if we choose strategy Alpha is:
                # P(success | i, j, Alpha) = P(success | i-1, j+1) * 0.6 + P(success | i-1, j-1) * 0.4
                # Each number from the problem description is used here:
                # Fee = £1 (so capital becomes j-1 before return)
                # Return on success = £2 (so final capital is j-1+2 = j+1)
                # Success probability = 60% (0.6)
                # Return on failure = £0 (so final capital is j-1+0 = j-1)
                # Failure probability = 40% (0.4)
                
                # Check bounds to avoid index out of range
                p_success_outcome = dp[i-1][j+1] if (j+1) <= MAX_CAPITAL else 0.0
                p_fail_outcome = dp[i-1][j-1] if (j-1) >= 0 else 0.0
                prob_alpha = 0.6 * p_success_outcome + 0.4 * p_fail_outcome

            # --- Calculate probability of success if we choose Strategy Beta ---
            prob_beta = 0.0
            # We must be able to afford the £3 fee.
            if j >= 3:
                # The final equation for the probability of success from this state (i, j)
                # if we choose strategy Beta is:
                # P(success | i, j, Beta) = P(success | i-1, j+12) * 0.2 + P(success | i-1, j-3) * 0.8
                # Each number from the problem description is used here:
                # Fee = £3 (so capital becomes j-3 before return)
                # Return on success = £15 (so final capital is j-3+15 = j+12)
                # Success probability = 20% (0.2)
                # Return on failure = £0 (so final capital is j-3+0 = j-3)
                # Failure probability = 80% (0.8)

                # Check bounds to avoid index out of range
                p_success_outcome = dp[i-1][j+12] if (j+12) <= MAX_CAPITAL else 0.0
                p_fail_outcome = dp[i-1][j-3] if (j-3) >= 0 else 0.0
                prob_beta = 0.2 * p_success_outcome + 0.8 * p_fail_outcome
            
            # The optimal strategy is to choose the action that maximizes the probability.
            dp[i][j] = max(prob_alpha, prob_beta)

    # The final answer is the probability of success starting with M capital and n trades.
    return dp[n][M]

if __name__ == '__main__':
    # Example usage of the function.
    # Let's use M=20, n=10. We want to reach £40.
    initial_investment = 20
    num_trades = 10
    
    # Run the algorithm
    result_probability = solve_trading_problem(initial_investment, num_trades)
    
    # Print the result
    print(f"Initial Investment (M): £{initial_investment}")
    print(f"Number of Trades (n): {num_trades}")
    print(f"Target Investment (2M): £{2 * initial_investment}")
    print(f"Maximum Probability of Success: {result_probability:.10f}")

    # Another interesting case: M=10, n=10. Target £20.
    # To get a net profit of £10 in 10 trades, one possible way is 10 successful Alpha trades.
    # The probability of that specific path is 0.6^10. Let's see what the DP finds.
    initial_investment_2 = 10
    num_trades_2 = 10
    result_probability_2 = solve_trading_problem(initial_investment_2, num_trades_2)
    print("\n--- Another Example ---")
    print(f"Initial Investment (M): £{initial_investment_2}")
    print(f"Number of Trades (n): {num_trades_2}")
    print(f"Target Investment (2M): £{2 * initial_investment_2}")
    print(f"Maximum Probability of Success: {result_probability_2:.10f}")
    print(f"(For comparison, probability of 10 straight Alpha successes is: {0.6**10:.10f})")

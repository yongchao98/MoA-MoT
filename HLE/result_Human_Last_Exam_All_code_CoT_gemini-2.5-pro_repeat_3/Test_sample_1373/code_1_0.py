def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    to exactly 2M in n trades using dynamic programming.
    """
    TARGET_MONEY = 2 * M
    
    # The maximum possible money you can have is M + n * 12 (n successful Beta trades).
    # Any state beyond this is unreachable, but we need to account for it in our table size.
    MAX_MONEY = M + 12 * n
    
    # dp[l][m] = max probability of success, given l trades are done and we have m money.
    # Initialize a 2D list for the DP table with all probabilities as 0.0.
    dp = [[0.0 for _ in range(MAX_MONEY + 1)] for _ in range(n + 1)]

    # Base case: at the end of n trades (l=n).
    # If we have exactly TARGET_MONEY, the probability of success is 1.
    if TARGET_MONEY <= MAX_MONEY:
        dp[n][TARGET_MONEY] = 1.0

    # Iterate backwards from l = n-1 down to 0 trades.
    for l in range(n - 1, -1, -1):
        # Iterate through all possible money amounts m for the current number of trades l.
        # The loop must go up to MAX_MONEY because a state at l can transition from a
        # state with more money at l+1 (e.g., m-3 at l+1 can lead to m at l).
        for m in range(MAX_MONEY + 1):
            
            # --- Calculate probability of success from this state (m, l) for each strategy ---

            # Strategy Alpha: Fee £1, Win +£2 (net +£1), Loss +£0 (net -£1)
            p_alpha = 0.0
            if m >= 1:
                # Get probability from the next state (l+1) after a win
                prob_alpha_win = dp[l+1][m+1] if m + 1 <= MAX_MONEY else 0.0
                # Get probability from the next state (l+1) after a loss
                prob_alpha_lose = dp[l+1][m-1]
                # The equation for Alpha's success probability
                p_alpha = 0.6 * prob_alpha_win + 0.4 * prob_alpha_lose

            # Strategy Beta: Fee £3, Win +£15 (net +£12), Loss +£0 (net -£3)
            p_beta = 0.0
            if m >= 3:
                # Get probability from the next state (l+1) after a win
                prob_beta_win = dp[l+1][m+12] if m + 12 <= MAX_MONEY else 0.0
                # Get probability from the next state (l+1) after a loss
                prob_beta_lose = dp[l+1][m-3]
                # The equation for Beta's success probability
                p_beta = 0.2 * prob_beta_win + 0.8 * prob_beta_lose

            # --- Determine the optimal strategy and store its probability in dp[l][m] ---
            # We must make a trade, so we choose the best option available.
            
            if m < 1:
                # Can't afford any trade, so probability of success from here is 0.
                dp[l][m] = 0.0
            elif m < 3:
                # Can only afford Alpha.
                dp[l][m] = p_alpha
            else: # m >= 3
                # Can afford both, so choose the one with the maximum probability.
                dp[l][m] = max(p_alpha, p_beta)
                
    # The final answer is the probability at the start: 0 trades done, M money.
    final_probability = dp[0][M]
    return final_probability

if __name__ == '__main__':
    # --- Example values for M and n ---
    # You can change these values to test different scenarios.
    initial_investment_M = 10
    num_trades_n = 10

    # Calculate the probability
    result = solve_trading_problem(initial_investment_M, num_trades_n)

    # --- Print the final results ---
    print(f"Initial Investment (M): £{initial_investment_M}")
    print(f"Number of Trades (n): {num_trades_n}")
    print(f"Target Investment (2M): £{2 * initial_investment_M}")
    print(f"\nThe maximum probability of reaching the target is: {result:.6f}")

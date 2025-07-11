def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling the initial investment M
    in exactly n trades using dynamic programming.
    """

    # The maximum possible capital is the initial capital plus n times the max profit from one trade (£12).
    max_capital = M + 12 * n

    # dp[i][j]: max probability of success with i trades LEFT and j capital.
    # Dimensions: (n+1) for trades, (max_capital+1) for money.
    dp = [[0.0] * (max_capital + 1) for _ in range(n + 1)]

    # Base Case: 0 trades left.
    # Success is true only if capital is exactly 2*M.
    target_capital = 2 * M
    if target_capital <= max_capital:
        dp[0][target_capital] = 1.0

    # Fill DP table using the recurrence relation, from i=1 to n trades left.
    for i in range(1, n + 1):
        # Iterate over all possible capital amounts.
        for j in range(max_capital + 1):

            # --- Calculate probability of success from choosing Strategy Alpha ---
            prob_alpha = 0.0
            if j >= 1: # Must have £1 to invest.
                # Net capital change: +1 (success) or -1 (failure)
                capital_alpha_success = j + 1
                capital_alpha_fail = j - 1
                
                # Get probability from the next state (i-1 trades left), checking bounds.
                prob_from_alpha_success = dp[i-1][capital_alpha_success] if capital_alpha_success <= max_capital else 0.0
                prob_from_alpha_fail = dp[i-1][capital_alpha_fail] if capital_alpha_fail >= 0 else 0.0
                
                # Expected probability of success for Alpha
                prob_alpha = 0.6 * prob_from_alpha_success + 0.4 * prob_from_alpha_fail

            # --- Calculate probability of success from choosing Strategy Beta ---
            prob_beta = 0.0
            if j >= 3: # Must have £3 to invest.
                # Net capital change: +12 (success) or -3 (failure)
                capital_beta_success = j + 12
                capital_beta_fail = j - 3
                
                # Get probability from the next state (i-1 trades left), checking bounds.
                prob_from_beta_success = dp[i-1][capital_beta_success] if capital_beta_success <= max_capital else 0.0
                prob_from_beta_fail = dp[i-1][capital_beta_fail] if capital_beta_fail >= 0 else 0.0

                # Expected probability of success for Beta
                prob_beta = 0.2 * prob_from_beta_success + 0.8 * prob_from_beta_fail

            # --- Determine the optimal strategy for state (i, j) ---
            # The choice depends on which strategies are affordable. A trade must be made.
            if j >= 3:
                dp[i][j] = max(prob_alpha, prob_beta)
            elif j >= 1:
                dp[i][j] = prob_alpha
            else:
                # Cannot afford any trade, so probability of success from this state is 0.
                dp[i][j] = 0.0

    final_prob = dp[n][M]
    return final_prob

def main():
    """
    Main function to run the simulation with example values and print the output.
    """
    # Example values for M and n.
    # Note: For many M, n combinations, the probability might be 0.
    # M=19, n=3 is an example that yields a non-zero probability.
    M = 19
    n = 3
    
    # The subproblem T[m, l] corresponds to dp[n-l][m] in our backward formulation.
    # The recurrence relation for T[m, l] (prob of success with m capital and l trades executed) is:
    # T[m, l] depends on T[m', l+1], working backwards from l=n.
    # For a state (m, l) (money m, l trades done), the probability of final success is T[l][m].
    # This is calculated from states at trade l+1. Let's express this using T.
    
    print("The DP recurrence relation can be expressed as follows:")
    print("Let T(m, l) be the max probability of success with £m and l trades remaining.")
    print("\nThe recurrence is:")
    print("For 1 <= m < 3: T(m, l) = 0.60 * T(m + 1, l - 1) + 0.40 * T(m - 1, l - 1)")
    print("For m >= 3:      T(m, l) = max( (0.60*T(m+1, l-1) + 0.40*T(m-1, l-1)), (0.20*T(m+12, l-1) + 0.80*T(m-3, l-1)) )")
    print("\nBase Case: T(m, 0) = 1.0 if m = 2*M, and 0.0 otherwise.")
    
    print("\n" + "="*50)
    print(f"Calculating for M = {M} and n = {n} trades...")
    
    result = solve_trading_probability(M, n)
    
    print(f"\nThe initial capital is £{M}. The target capital is £{2*M}.")
    print(f"The maximum probability of reaching exactly £{2*M} in {n} trades is: {result}")
    print("="*50)

if __name__ == "__main__":
    main()
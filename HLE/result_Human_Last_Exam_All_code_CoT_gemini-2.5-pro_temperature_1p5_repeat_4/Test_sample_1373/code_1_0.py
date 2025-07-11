import collections

def solve_trading_probability():
    """
    Calculates the maximum probability of doubling an initial investment
    using dynamic programming.
    """
    # --- Problem Parameters ---
    M = 25  # Initial investment in GBP
    n = 5   # Exact number of trades to execute

    # The target is to double the initial investment
    target_m = 2 * M

    # --- Dynamic Programming Setup ---
    # dp_tables[l] will be a dictionary holding the probabilities for states
    # with l trades remaining. Using defaultdict is efficient for sparse states.
    dp_tables = [collections.defaultdict(float) for _ in range(n + 1)]

    # Base Case: l = 0 trades remaining.
    # Probability of success is 1.0 only if we have exactly the target amount.
    dp_tables[0][target_m] = 1.0

    # --- DP Calculation Loop ---
    # Iterate from l = 1 to n trades remaining.
    for l in range(1, n + 1):
        # Determine the range of money 'm' to calculate for.
        # An optimistic bound for the money we can have after (n-l) trades is M + 12*(n-l)
        # A pessimistic bound is M - 3*(n-l)
        # We'll iterate over a safe range around the target.
        # Minimum possible money to reach target: target_m - 12*l
        # Maximum possible money to reach target: target_m + 3*l
        min_m_range = max(0, target_m - 12 * l)
        max_m_range = target_m + 3 * l

        for m in range(min_m_range, max_m_range + 15): # Add buffer for m+12 lookups
            # Probability of success if we choose Strategy Alpha
            prob_alpha = 0.0
            if m >= 1:
                prob_alpha = (0.60 * dp_tables[l - 1][m + 1] +
                              0.40 * dp_tables[l - 1][m - 1])

            # Probability of success if we choose Strategy Beta
            prob_beta = 0.0
            if m >= 3:
                prob_beta = (0.20 * dp_tables[l - 1][m + 12] +
                             0.80 * dp_tables[l - 1][m - 3])

            # Store the maximum probability from the two choices
            dp_tables[l][m] = max(prob_alpha, prob_beta)

    # --- Final Result and Equation ---
    final_probability = dp_tables[n][M]

    print(f"Initial Investment (M): £{M}")
    print(f"Number of Trades (n): {n}")
    print(f"Target Investment (2M): £{target_m}\n")
    print(f"Maximum probability of success: {final_probability:.6f}\n")
    
    # Show the calculation for the final state (the first trade)
    print("--- Calculation for the First Trade Decision ---")
    
    # Get the probabilities from the two strategies for the initial state (n, M)
    prob_alpha_final = 0.0
    if M >= 1:
        prob_alpha_final = (0.60 * dp_tables[n-1][M+1] + 
                            0.40 * dp_tables[n-1][M-1])
                            
    prob_beta_final = 0.0
    if M >= 3:
        prob_beta_final = (0.20 * dp_tables[n-1][M+12] + 
                           0.80 * dp_tables[n-1][M-3])

    if prob_alpha_final >= prob_beta_final:
        print("Optimal initial choice: Strategy Alpha")
        print("Equation: P(start) = 0.60 * P(state after win) + 0.40 * P(state after loss)")
        print(f"Breakdown: {final_probability:.6f} = 0.60 * {dp_tables[n-1][M+1]:.6f} + 0.40 * {dp_tables[n-1][M-1]:.6f}")
    else:
        print("Optimal initial choice: Strategy Beta")
        print("Equation: P(start) = 0.20 * P(state after win) + 0.80 * P(state after loss)")
        print(f"Breakdown: {final_probability:.6f} = 0.20 * {dp_tables[n-1][M+12]:.6f} + 0.80 * {dp_tables[n-1][M-3]:.6f}")

# Execute the algorithm
solve_trading_probability()

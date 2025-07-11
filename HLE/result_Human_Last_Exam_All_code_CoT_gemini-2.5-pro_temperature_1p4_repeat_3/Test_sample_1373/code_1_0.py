def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    Let T[l][m] be the maximum probability of success (ending with 2M after n trades)
    given we have completed l trades and have m pounds.
    """
    
    # Determine the maximum money we might need to track.
    # The highest possible amount is starting with M and always getting the
    # successful outcome of Strategy Beta (+£12 profit per trade).
    max_m_attainable = M + 12 * n
    # The target money is 2*M.
    max_m_target = 2 * M
    
    # max_money determines the size of our DP table's money dimension.
    # We add a buffer because the largest state we might query is m+12.
    max_money = max(max_m_attainable, max_m_target) + 13

    # dp[l][m] corresponds to the problem's T[m, l].
    dp = [[0.0] * max_money for _ in range(n + 1)]

    # Base Case: At the end of n trades (l=n).
    # Success is only achieved if we have exactly 2*M pounds.
    if 2 * M < max_money:
        dp[n][2 * M] = 1.0

    # Fill the DP table by iterating backwards from trade l = n-1 to 0.
    for l in range(n - 1, -1, -1):
        # Iterate over all possible money amounts at this stage.
        for m in range(max_money):
            prob_alpha = 0.0
            prob_beta = 0.0

            # --- Strategy Alpha (Cost: £1) ---
            if m >= 1:
                # To calculate T[m,l], we look at the probabilities at trade l+1.
                # Success (60%): money becomes m - 1 + 2 = m + 1. Prob from here is dp[l+1][m+1].
                # Failure (40%): money becomes m - 1 + 0 = m - 1. Prob from here is dp[l+1][m-1].
                prob_alpha = 0.60 * dp[l + 1][m + 1] + 0.40 * dp[l + 1][m - 1]

            # --- Strategy Beta (Cost: £3) ---
            if m >= 3:
                # Success (20%): money becomes m - 3 + 15 = m + 12. Prob is dp[l+1][m+12].
                # Failure (80%): money becomes m - 3 + 0 = m - 3. Prob is dp[l+1][m-3].
                prob_beta = 0.20 * dp[l + 1][m + 12] + 0.80 * dp[l + 1][m - 3]

            # The optimal strategy at state (l, m) is to choose the action
            # that maximizes the probability of future success.
            if m < 1:
                dp[l][m] = 0.0  # Cannot afford any trade.
            elif m < 3:
                dp[l][m] = prob_alpha  # Can only afford Alpha.
            else:
                dp[l][m] = max(prob_alpha, prob_beta)  # Can afford both, choose the best.
    
    # --- Final Result and Explanation ---
    final_prob = dp[0][M]
    print(f"Goal: Double initial investment from £{M} to £{2*M} in exactly {n} trades.\n")
    
    print(f"The calculation for the initial state T[m={M}, l=0]:")

    if M >= 1:
        prob_alpha_final = 0.60 * dp[1][M + 1] + 0.40 * dp[1][M - 1]
        print(f"  Prob. from Alpha = 0.60 * T[m={M+1}, l=1] + 0.40 * T[m={M-1}, l=1]")
        print(f"                   = 0.60 * {dp[1][M + 1]:.6f} + 0.40 * {dp[1][M - 1]:.6f}")
        print(f"                   = {prob_alpha_final:.6f}")
    else:
        print("  Strategy Alpha is not affordable at the start.")

    if M >= 3:
        prob_beta_final = 0.20 * dp[1][M + 12] + 0.80 * dp[1][M - 3]
        print(f"  Prob. from Beta  = 0.20 * T[m={M+12}, l=1] + 0.80 * T[m={M-3}, l=1]")
        print(f"                   = 0.20 * {dp[1][M + 12]:.6f} + 0.80 * {dp[1][M - 3]:.6f}")
        print(f"                   = {prob_beta_final:.6f}")
    else:
        print("  Strategy Beta is not affordable at the start.")

    print("\nOptimal Choice at Start:")
    if M >= 3:
        if prob_alpha_final > prob_beta_final:
            print("  The optimal first move is Strategy Alpha.")
        elif prob_beta_final > prob_alpha_final:
             print("  The optimal first move is Strategy Beta.")
        else:
             print("  Both strategies yield the same probability; either can be chosen.")
    elif M >= 1:
        print("  The only available move is Strategy Alpha.")
    else:
        print("  No trades are possible.")


    print("\n------------------------------------------------------------")
    print(f"The maximum probability of success is: {final_prob:.6f}")
    print("------------------------------------------------------------")


# --- Parameters for the simulation ---
# Initial Money in GBP
initial_money = 25
# Number of Trades to be executed
num_trades = 10

solve_trading_problem(initial_money, num_trades)
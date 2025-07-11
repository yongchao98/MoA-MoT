def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.
    
    The subproblem is defined as:
    dp[l][m] = Maximum probability of success with 'l' trades remaining and 'm' money.
    """
    
    # The target amount of money to reach is 2*M
    target_money = 2 * M

    # The amount of money can increase. The maximum increase per trade is +£12 (from a successful Beta trade).
    # The maximum possible money after n trades is M + n * 12.
    # The DP table's money dimension must be large enough to hold the target and all possible intermediate values.
    # A safe upper limit for the money dimension is max(target_money, M + n * 12) plus a buffer for indexing (e.g., m+12).
    money_limit = max(target_money, M + n * 12) + 15
    
    # Initialize the DP table.
    # Dimensions: (n + 1) rows for trades remaining, (money_limit) columns for money.
    dp = [[0.0 for _ in range(money_limit)] for _ in range(n + 1)]

    # Base Case: l = 0 (no trades left)
    # If we have exactly the target amount, the probability of success is 1. Otherwise, it's 0.
    if target_money < money_limit:
        dp[0][target_money] = 1.0

    # Fill the DP table bottom-up using the recurrence relation.
    # Iterate from l = 1 to n trades remaining.
    for l in range(1, n + 1):
        # Iterate through all relevant money amounts.
        for m in range(money_limit):
            
            # --- Calculate probability if we choose Strategy Alpha ---
            prob_alpha = 0.0
            if m >= 1: # Must have at least £1 to pay the fee.
                # Future states after an Alpha trade: m+1 (success) or m-1 (failure).
                # We look up the probabilities from the already computed l-1 subproblems.
                prob_success_alpha = dp[l-1][m + 1] if (m + 1) < money_limit else 0.0
                prob_fail_alpha = dp[l-1][m - 1] if (m - 1) >= 0 else 0.0
                prob_alpha = 0.60 * prob_success_alpha + 0.40 * prob_fail_alpha
            
            # --- Calculate probability if we choose Strategy Beta ---
            prob_beta = 0.0
            if m >= 3: # Must have at least £3 to pay the fee.
                # Future states after a Beta trade: m+12 (success) or m-3 (failure).
                prob_success_beta = dp[l-1][m + 12] if (m + 12) < money_limit else 0.0
                prob_fail_beta = dp[l-1][m - 3] if (m - 3) >= 0 else 0.0
                prob_beta = 0.20 * prob_success_beta + 0.80 * prob_fail_beta
            
            # The optimal strategy is to choose the action with the maximum probability of success.
            dp[l][m] = max(prob_alpha, prob_beta)
            
    # The final answer is the probability starting with M money and n trades remaining.
    final_prob = dp[n][M]
    
    # As requested, print the breakdown of the final calculation.
    print(f"Calculating the optimal probability for an initial investment of £{M} and {n} trades:")
    print("-" * 60)
    print(f"The goal is to reach £{target_money}.")
    print(f"The desired value is P(trades_left={n}, money={M}).\n")

    # Retrieve values needed to show the final calculation.
    prob_alpha_final, prob_beta_final = 0.0, 0.0
    
    if M >= 1:
        val_alpha_succ = dp[n-1][M+1]
        val_alpha_fail = dp[n-1][M-1]
        prob_alpha_final = 0.60 * val_alpha_succ + 0.40 * val_alpha_fail
        print("Final choice - Strategy Alpha:")
        print(f"  P_alpha = 0.60 * P(trades_left={n-1}, money={M+1}) + 0.40 * P(trades_left={n-1}, money={M-1})")
        print(f"  P_alpha = 0.60 * {val_alpha_succ:.4f} + 0.40 * {val_alpha_fail:.4f}")
        print(f"  P_alpha = {prob_alpha_final:.4f}\n")
    else:
        print("Final choice - Strategy Alpha not possible (not enough money).\n")

    if M >= 3:
        val_beta_succ = dp[n-1][M+12]
        val_beta_fail = dp[n-1][M-3]
        prob_beta_final = 0.20 * val_beta_succ + 0.80 * val_beta_fail
        print("Final choice - Strategy Beta:")
        print(f"  P_beta = 0.20 * P(trades_left={n-1}, money={M+12}) + 0.80 * P(trades_left={n-1}, money={M-3})")
        print(f"  P_beta = 0.20 * {val_beta_succ:.4f} + 0.80 * {val_beta_fail:.4f}")
        print(f"  P_beta = {prob_beta_final:.4f}\n")
    else:
        print("Final choice - Strategy Beta not possible (not enough money).\n")

    print("-" * 60)
    print(f"The optimal strategy chooses the maximum of these probabilities.")
    print(f"Final Optimal Probability = max({prob_alpha_final:.4f}, {prob_beta_final:.4f})")
    print(f"Result: {final_prob:.4f}")

# --- Main execution block ---
if __name__ == '__main__':
    # You can change these values to test different scenarios.
    # For M=2, n=2, the only path to the target of £4 is two successful Alpha trades (0.6*0.6=0.36).
    initial_investment_M = 2
    number_of_trades_n = 2
    
    solve_trading_problem(initial_investment_M, number_of_trades_n)
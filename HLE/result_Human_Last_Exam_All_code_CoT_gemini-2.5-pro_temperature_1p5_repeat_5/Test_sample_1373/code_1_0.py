def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.
    """
    if M <= 0 or n < 0:
        print("Initial investment and number of trades must be non-negative.")
        return

    # Maximum possible money after n trades with maximum profit per trade (+12)
    max_money = M + 12 * n
    target_money = 2 * M

    # dp[k][m] = max probability of success with k trades remaining and m pounds.
    # Dimensions: (n+1) trades x (max_money+1) money amounts
    dp = [[0.0 for _ in range(max_money + 1)] for _ in range(n + 1)]

    # Base case: k = 0 trades remaining.
    # Success is 100% (1.0) if money is exactly target_money, 0% otherwise.
    if target_money <= max_money:
        dp[0][target_money] = 1.0

    # Fill the DP table bottom-up (from k=1 to n trades remaining)
    for k in range(1, n + 1):
        for m in range(max_money + 1):
            prob_alpha = 0.0
            prob_beta = 0.0

            # Option 1: Strategy Alpha (Cost: £1)
            if m >= 1:
                # Value from success state (m+1)
                p_succ_alpha = dp[k - 1][m + 1] if m + 1 <= max_money else 0.0
                # Value from failure state (m-1)
                p_fail_alpha = dp[k - 1][m - 1] if m - 1 >= 0 else 0.0
                prob_alpha = 0.6 * p_succ_alpha + 0.4 * p_fail_alpha

            # Option 2: Strategy Beta (Cost: £3)
            if m >= 3:
                # Value from success state (m+12)
                p_succ_beta = dp[k - 1][m + 12] if m + 12 <= max_money else 0.0
                # Value from failure state (m-3)
                p_fail_beta = dp[k - 1][m - 3] if m - 3 >= 0 else 0.0
                prob_beta = 0.2 * p_succ_beta + 0.8 * p_fail_beta
            
            # The optimal strategy is the one with the maximum probability of success
            dp[k][m] = max(prob_alpha, prob_beta)

    final_probability = dp[n][M]
    
    # Output the details of the calculation for the first trade (at state n, M)
    print(f"To find the optimal probability starting with M={M} GBP and n={n} trades:")
    print("-" * 30)

    can_use_alpha = M >= 1
    can_use_beta = M >= 3

    prob_alpha_final = 0.0
    if can_use_alpha:
        p_succ_alpha = dp[n - 1][M + 1] if M + 1 <= max_money else 0
        p_fail_alpha = dp[n - 1][M - 1] if M - 1 >= 0 else 0
        prob_alpha_final = 0.6 * p_succ_alpha + 0.4 * p_fail_alpha
        print("Strategy Alpha Calculation for the first trade:")
        print(f"Prob = 0.6 * P(success | {M+1} money, {n-1} trades left) + 0.4 * P(success | {M-1} money, {n-1} trades left)")
        print(f"     = 0.6 * {p_succ_alpha:.4f} + 0.4 * {p_fail_alpha:.4f}")
        print(f"     = {prob_alpha_final:.4f}\n")

    prob_beta_final = 0.0
    if can_use_beta:
        p_succ_beta = dp[n-1][M+12] if M+12 <= max_money else 0
        p_fail_beta = dp[n-1][M-3] if M-3 >= 0 else 0
        prob_beta_final = 0.2 * p_succ_beta + 0.8 * p_fail_beta
        print("Strategy Beta Calculation for the first trade:")
        print(f"Prob = 0.2 * P(success | {M+12} money, {n-1} trades left) + 0.8 * P(success | {M-3} money, {n-1} trades left)")
        print(f"     = 0.2 * {p_succ_beta:.4f} + 0.8 * {p_fail_beta:.4f}")
        print(f"     = {prob_beta_final:.4f}\n")

    if not can_use_alpha and not can_use_beta:
        print("No trades are possible with the initial capital.")
    elif prob_alpha_final > prob_beta_final:
        print("Optimal first move is Strategy Alpha.")
    elif prob_beta_final > prob_alpha_final:
        print("Optimal first move is Strategy Beta.")
    else:
        if can_use_alpha or can_use_beta:
            print("Both strategies yield the same probability, so either can be chosen.")
        
    print("-" * 30)
    print(f"The maximum probability of reaching exactly {target_money} GBP is: {final_probability}")

# --- You can change these input values to test the algorithm ---
initial_investment_M = 20
number_of_trades_n = 5
solve_trading_problem(initial_investment_M, number_of_trades_n)
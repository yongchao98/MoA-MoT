def solve_trading_probability():
    """
    Calculates the maximum probability of doubling an initial investment M
    by making exactly n trades.
    """
    # Problem parameters (example values, as none were provided)
    M = 20  # Initial investment
    n = 10   # Number of trades

    # The maximum possible amount of money we can have.
    # Start with M, max profit per trade is 12 (from Beta's 15 profit - 3 fee).
    MAX_MONEY = M + n * 12

    # DP table: T[k][m] = max probability of success with k trades remaining and m money.
    # Dimensions: (n+1) rows for trades, (MAX_MONEY+1) cols for money.
    T = [[0.0 for _ in range(MAX_MONEY + 1)] for _ in range(n + 1)]

    # Base case: k = 0 trades remaining
    # If we have exactly 2M pounds, success probability is 1.0.
    target_money = 2 * M
    if target_money <= MAX_MONEY:
        T[0][target_money] = 1.0

    # Fill the DP table using the recurrence relation
    # k = number of trades remaining
    for k in range(1, n + 1):
        # m = current amount of money
        for m in range(MAX_MONEY + 1):
            
            # Probability from Strategy Alpha (if affordable)
            prob_alpha = 0.0
            if m >= 1:
                # Net profit is 1, so money becomes m+1
                p_success_alpha = T[k - 1][m + 1] if (m + 1) <= MAX_MONEY else 0.0
                # Net loss is 1, so money becomes m-1
                p_failure_alpha = T[k - 1][m - 1] if (m - 1) >= 0 else 0.0
                prob_alpha = 0.60 * p_success_alpha + 0.40 * p_failure_alpha

            # Probability from Strategy Beta (if affordable)
            prob_beta = 0.0
            if m >= 3:
                # Net profit is 12, so money becomes m+12
                p_success_beta = T[k - 1][m + 12] if (m + 12) <= MAX_MONEY else 0.0
                # Net loss is 3, so money becomes m-3
                p_failure_beta = T[k - 1][m - 3] if (m - 3) >= 0 else 0.0
                prob_beta = 0.20 * p_success_beta + 0.80 * p_failure_beta

            # The optimal strategy is the one with the maximum probability of success
            optimal_prob = 0.0
            if m >= 3:
                optimal_prob = max(prob_alpha, prob_beta)
            elif m >= 1:
                optimal_prob = prob_alpha
            
            T[k][m] = optimal_prob
            
    # The final answer is T(n, M)
    final_probability = T[n][M]
    
    print("This program calculates the probability of success using dynamic programming.")
    print(f"Initial investment M = {M} GBP, Number of trades n = {n}.")
    print(f"Target investment 2M = {2*M} GBP.")
    print("\nLet T(k, m) be the max probability of success with k trades remaining and money m.")
    print("The recurrence relation is T(k, m) = max(Prob_Alpha, Prob_Beta) where:")
    
    print("Prob_Alpha = 0.60 * T(k-1, m - 1 + 2) + 0.40 * T(k-1, m - 1 + 0)")
    print("Prob_Beta  = 0.20 * T(k-1, m - 3 + 15) + 0.80 * T(k-1, m - 3 + 0)")
    print("\nSimplified:")
    print("Prob_Alpha = 0.60 * T(k-1, m + 1) + 0.40 * T(k-1, m - 1)       (Requires m >= 1)")
    print("Prob_Beta  = 0.20 * T(k-1, m + 12) + 0.80 * T(k-1, m - 3)      (Requires m >= 3)")
          
    print(f"\nThe optimal probability of reaching exactly {2*M} GBP is: {final_probability}")

if __name__ == '__main__':
    solve_trading_probability()
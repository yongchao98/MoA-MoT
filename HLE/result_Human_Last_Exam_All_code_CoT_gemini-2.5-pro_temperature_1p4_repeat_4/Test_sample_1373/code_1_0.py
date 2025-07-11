def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    The subproblem T[m, l] (or dp[l][m] in code) represents the maximum probability 
    of success given that we start trade 'l' with 'm' GBP. The algorithm works 
    backward from the final trade 'n'.

    Args:
        M (int): Initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        float: The maximum probability of reaching 2*M GBP.
    """

    # --- Strategy Parameters ---
    # Alpha: cost 1, 60% for +1 net, 40% for -1 net
    alpha_cost = 1
    alpha_prob_success = 0.6
    alpha_net_gain = 1  # (2 return - 1 fee)
    alpha_net_loss = -1 # (0 return - 1 fee)

    # Beta: cost 3, 20% for +12 net, 80% for -3 net
    beta_cost = 3
    beta_prob_success = 0.2
    beta_net_gain = 12 # (15 return - 3 fee)
    beta_net_loss = -3 # (0 return - 3 fee)

    # --- DP Table Initialization ---
    # The maximum possible capital after n trades with max gain of +12.
    # We add a buffer (equal to the max gain) to avoid index out of bounds
    # when calculating future states like m + 12.
    max_capital = M + beta_net_gain * n + beta_net_gain
    target_capital = 2 * M

    # dp[l][m]: max probability of success starting trade l with capital m
    # Dimensions: (n+1) trades x (max_capital+1) capital amounts
    dp = [[0.0 for _ in range(max_capital + 1)] for _ in range(n + 1)]

    # --- Base Case ---
    # At the end of n trades (l=n), success is 1.0 only if capital is exactly 2*M.
    if target_capital >= 0 and target_capital <= max_capital:
        dp[n][target_capital] = 1.0

    # --- DP Calculation (Iterating Backwards) ---
    for l in range(n - 1, -1, -1):  # Iterate from trade n-1 down to 0
        for m in range(max_capital + 1): # Iterate through all possible capital amounts
            
            # Probability of success if we choose Alpha
            prob_alpha = 0.0
            if m >= alpha_cost:
                # State after successful Alpha trade
                m_succ_alpha = m + alpha_net_gain
                prob_from_succ_alpha = dp[l+1][m_succ_alpha]

                # State after failed Alpha trade
                m_fail_alpha = m + alpha_net_loss
                prob_from_fail_alpha = 0.0
                if m_fail_alpha >= 0:
                    prob_from_fail_alpha = dp[l+1][m_fail_alpha]
                
                prob_alpha = (alpha_prob_success * prob_from_succ_alpha) + \
                             ((1 - alpha_prob_success) * prob_from_fail_alpha)

            # Probability of success if we choose Beta
            prob_beta = 0.0
            if m >= beta_cost:
                # State after successful Beta trade
                m_succ_beta = m + beta_net_gain
                prob_from_succ_beta = dp[l+1][m_succ_beta]

                # State after failed Beta trade
                m_fail_beta = m + beta_net_loss
                prob_from_fail_beta = 0.0
                if m_fail_beta >= 0:
                    prob_from_fail_beta = dp[l+1][m_fail_beta]

                prob_beta = (beta_prob_success * prob_from_succ_beta) + \
                            ((1 - beta_prob_success) * prob_from_fail_beta)
            
            # The optimal strategy is to choose the one with the max probability
            dp[l][m] = max(prob_alpha, prob_beta)

    # The final answer is the probability at the start (trade 0, capital M)
    final_probability = dp[0][M]
    return final_probability

# --- Example Usage ---
initial_investment_M = 25
num_trades_n = 5

# Calculate and print the result
result = solve_trading_problem(initial_investment_M, num_trades_n)
print(f"Starting with M = {initial_investment_M} GBP and executing n = {num_trades_n} trades.")
print(f"The maximum probability of doubling the investment to {2 * initial_investment_M} GBP is: {result}")
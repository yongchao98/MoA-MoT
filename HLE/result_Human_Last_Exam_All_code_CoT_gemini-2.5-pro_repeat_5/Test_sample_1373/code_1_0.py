import sys

# Increase recursion limit for deep traces, although our DP is iterative.
sys.setrecursionlimit(2000)

def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.
    """
    if M <= 0 or n <= 0:
        print("Initial investment and number of trades must be positive.")
        return

    # The maximum possible amount of money can be estimated.
    # If we always use Beta and win, money increases by 12 each time.
    # Max money is roughly M + 12*n. We add a buffer for safety in indexing.
    # The max index accessed is m+12, so the size needs to be M + 12*n + 13.
    max_money = M + 12 * n + 13

    # dp[l][m] = max probability of success with l trades remaining and m money.
    # We use a list of dictionaries for sparse storage, which can be more memory-efficient
    # if the number of reachable money states is small. A 2D list is also fine.
    dp = [{} for _ in range(n + 1)]

    # Base Case: l = 0 trades remaining.
    # Success only if money is exactly 2*M.
    target_money = 2 * M
    if target_money < max_money:
        dp[0][target_money] = 1.0

    # Fill DP table bottom-up, from l=1 to l=n trades remaining.
    for l in range(1, n + 1):
        # The range of money we can have after (n-l) trades is bounded.
        # Min: M - 3*(n-l), Max: M + 12*(n-l). We iterate over a safe superset.
        # This optimization is not strictly needed for correctness but improves performance.
        min_reachable_m = max(0, M - 3 * l)
        max_reachable_m = M + 12 * l
        
        for m in range(min_reachable_m, max_reachable_m + 1):
            if m >= max_money: continue

            # Probability of success if we choose Strategy Alpha
            p_alpha = 0.0
            if m >= 1:
                prob_succ_alpha = dp[l - 1].get(m + 1, 0.0)
                prob_fail_alpha = dp[l - 1].get(m - 1, 0.0)
                p_alpha = 0.6 * prob_succ_alpha + 0.4 * prob_fail_alpha

            # Probability of success if we choose Strategy Beta
            p_beta = 0.0
            if m >= 3:
                prob_succ_beta = dp[l - 1].get(m + 12, 0.0)
                prob_fail_beta = dp[l - 1].get(m - 3, 0.0)
                p_beta = 0.2 * prob_succ_beta + 0.8 * prob_fail_beta
            
            # The optimal strategy maximizes the probability
            final_prob_for_state = max(p_alpha, p_beta)
            if final_prob_for_state > 0:
                dp[l][m] = final_prob_for_state

    # --- Final Result and Equation Printing ---
    final_answer = dp[n].get(M, 0.0)

    print(f"Initial investment M = £{M}, Number of trades n = {n}")
    print(f"Target to achieve: £{target_money}")
    print("-" * 50)
    print("Let T(m, l) be the max probability of success with £m and l trades remaining.")
    print(f"The goal is to compute T({M}, {n}).")
    print("-" * 50)

    # Re-calculate the top-level decision to show the final equation
    p_alpha_final = 0.0
    val_alpha_succ, val_alpha_fail = 0.0, 0.0
    if M >= 1:
        val_alpha_succ = dp[n - 1].get(M + 1, 0.0)
        val_alpha_fail = dp[n - 1].get(M - 1, 0.0)
        p_alpha_final = 0.6 * val_alpha_succ + 0.4 * val_alpha_fail

    p_beta_final = 0.0
    val_beta_succ, val_beta_fail = 0.0, 0.0
    if M >= 3:
        val_beta_succ = dp[n - 1].get(M + 12, 0.0)
        val_beta_fail = dp[n - 1].get(M - 3, 0.0)
        p_beta_final = 0.2 * val_beta_succ + 0.8 * val_beta_fail

    print(f"For the first trade (when l={n}, m={M}), we evaluate the two strategies:")

    print("\n1. Probability from Strategy Alpha:")
    if M >= 1:
        print(f"   = 0.6 * T({M+1}, {n-1}) + 0.4 * T({M-1}, {n-1})")
        print(f"   = 0.6 * {val_alpha_succ:.6f} + 0.4 * {val_alpha_fail:.6f}")
        print(f"   = {p_alpha_final:.6f}")
    else:
        print("   Not possible (not enough money for £1 fee).")

    print("\n2. Probability from Strategy Beta:")
    if M >= 3:
        print(f"   = 0.2 * T({M+12}, {n-1}) + 0.8 * T({M-3}, {n-1})")
        print(f"   = 0.2 * {val_beta_succ:.6f} + 0.8 * {val_beta_fail:.6f}")
        print(f"   = {p_beta_final:.6f}")
    else:
        print("   Not possible (not enough money for £3 fee).")

    print("-" * 50)
    
    choice = "None"
    if p_alpha_final > p_beta_final:
        choice = "Alpha"
    elif p_beta_final > p_alpha_final:
        choice = "Beta"
    elif p_alpha_final > 0: # If they are equal and non-zero
        choice = "Alpha or Beta"

    print(f"Optimal first move is Strategy: {choice}")
    print(f"\nThe maximum probability T({M}, {n}) of doubling the investment is: {final_answer:.7f}")


if __name__ == '__main__':
    # Example values for M and n
    initial_investment_M = 25
    number_of_trades_n = 10
    solve_trading_probability(initial_investment_M, number_of_trades_n)

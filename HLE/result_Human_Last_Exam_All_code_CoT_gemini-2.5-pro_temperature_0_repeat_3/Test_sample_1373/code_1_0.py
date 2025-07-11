import sys

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
        # If M is the target, probability is 1 if n=0, 0 otherwise.
        # But problem implies trading, so n>0.
        if M > 0 and 2*M == M: # M=0 case
             print("The final optimal probability of success is: 1.0")
        else:
             print("The final optimal probability of success is: 0.0")
        return

    # The maximum amount of money we need to track is 2*M.
    # Any amount > 2*M is a failure state (prob=0).
    max_money = 2 * M

    # We use a space-optimized DP approach.
    # 'prev_dp' stores the probabilities for the next trade step (l+1).
    # Initialize for the base case l = n.
    # At l=n, probability is 1.0 only if money is exactly 2*M.
    prev_dp = [0.0] * (max_money + 1)
    prev_dp[max_money] = 1.0

    # Iterate backwards from trade l = n-1 down to l = 1.
    # We stop at l=1 because we want to show the detailed calculation for l=0.
    for l in range(n - 1, 0, -1):
        current_dp = [0.0] * (max_money + 1)
        for m in range(max_money + 1):
            # Probability of success if we choose Strategy Alpha
            p_alpha = 0.0
            if m >= 1:
                prob_success = 0.60 * (prev_dp[m + 1] if m + 1 <= max_money else 0.0)
                prob_failure = 0.40 * (prev_dp[m - 1] if m - 1 >= 0 else 0.0)
                p_alpha = prob_success + prob_failure

            # Probability of success if we choose Strategy Beta
            p_beta = 0.0
            if m >= 3:
                prob_success = 0.20 * (prev_dp[m + 12] if m + 12 <= max_money else 0.0)
                prob_failure = 0.80 * (prev_dp[m - 3] if m - 3 >= 0 else 0.0)
                p_beta = prob_success + prob_failure
            
            # The optimal choice maximizes the probability
            current_dp[m] = max(p_alpha, p_beta)
        
        prev_dp = current_dp

    # Now, 'prev_dp' holds the probabilities T[m, 1].
    # We can calculate the final answer for T[M, 0] and show the equation.
    
    # --- Calculate P_alpha for the initial state (M, 0) ---
    p_alpha_final = 0.0
    # T[M+1, 1]
    prob_alpha_succ_state = prev_dp[M + 1] if M + 1 <= max_money else 0.0
    # T[M-1, 1]
    prob_alpha_fail_state = prev_dp[M - 1] if M - 1 >= 0 else 0.0
    if M >= 1:
        p_alpha_final = 0.60 * prob_alpha_succ_state + 0.40 * prob_alpha_fail_state

    # --- Calculate P_beta for the initial state (M, 0) ---
    p_beta_final = 0.0
    # T[M+12, 1]
    prob_beta_succ_state = prev_dp[M + 12] if M + 12 <= max_money else 0.0
    # T[M-3, 1]
    prob_beta_fail_state = prev_dp[M - 3] if M - 3 >= 0 else 0.0
    if M >= 3:
        p_beta_final = 0.20 * prob_beta_succ_state + 0.80 * prob_beta_fail_state

    # --- Output the final equation breakdown ---
    print(f"Calculating the optimal probability for initial state (Money={M}, Trades Left={n})...")
    print("This is T[M, 0] = max(P_alpha, P_beta), where P_alpha and P_beta are the probabilities of success after one trade.")
    print("\n" + "="*60)
    print("Final Equation Breakdown:")
    print("="*60)
    
    print(f"\n1. Strategy Alpha (requires £1):")
    print(f"   P_alpha = 0.60 * T(m={M+1}, l=1) + 0.40 * T(m={M-1}, l=1)")
    print(f"   P_alpha = 0.60 * {prob_alpha_succ_state:.6f} + 0.40 * {prob_alpha_fail_state:.6f}")
    print(f"   P_alpha = {p_alpha_final:.6f}")

    print(f"\n2. Strategy Beta (requires £3):")
    print(f"   P_beta  = 0.20 * T(m={M+12}, l=1) + 0.80 * T(m={M-3}, l=1)")
    print(f"   P_beta  = 0.20 * {prob_beta_succ_state:.6f} + 0.80 * {prob_beta_fail_state:.6f}")
    print(f"   P_beta  = {p_beta_final:.6f}")

    final_probability = max(p_alpha_final, p_beta_final)
    
    print("\n" + "="*60)
    print("Optimal Choice:")
    print(f"T(M={M}, l=0) = max({p_alpha_final:.6f}, {p_beta_final:.6f})")
    print(f"\nThe final optimal probability of success is: {final_probability:.6f}")


# --- Example Usage ---
# You can change these values to test different scenarios.
initial_investment_M = 20
number_of_trades_n = 10

solve_trading_probability(initial_investment_M, number_of_trades_n)
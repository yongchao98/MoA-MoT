import collections

def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    by making exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.
    """
    if M <= 0 or n <= 0:
        print("Initial investment and number of trades must be positive.")
        # If the goal is 2*M=0 and we start with 0, prob is 1 if n=0, 0 otherwise.
        # This case is simplified here.
        if 2*M == M:
            print("Final probability: 1.0000")
        else:
            print("Final probability: 0.0000")
        return

    target_money = 2 * M

    # dp[m] will store the max probability of success with the current number of trades
    # remaining and a capital of m. We use a defaultdict to handle the sparse nature
    # of reachable money states efficiently.
    # Initialize for the base case: l = 0 trades remaining.
    dp = collections.defaultdict(float)
    dp[target_money] = 1.0

    # This will store the state for n-1 trades, needed to print the final equation.
    dp_n_minus_1 = None

    # Iterate backwards from l = 1 to n trades remaining.
    for l in range(1, n + 1):
        # Before computing the probabilities for l=n, save the l=n-1 state.
        if l == n:
            dp_n_minus_1 = dp.copy()

        new_dp = collections.defaultdict(float)
        
        # Determine the relevant range of money 'm' to calculate for.
        # A state 'm' with 'l' trades remaining is reachable from M after (n-l) trades.
        trades_done = n - l
        min_m_range = M - 3 * trades_done
        max_m_range = M + 12 * trades_done
        
        for m in range(min_m_range, max_m_range + 1):
            # Probability of success if we choose Strategy Alpha
            prob_alpha = 0.0
            if m >= 1:
                prob_alpha = 0.6 * dp[m + 1] + 0.4 * dp[m - 1]

            # Probability of success if we choose Strategy Beta
            prob_beta = 0.0
            if m >= 3:
                prob_beta = 0.2 * dp[m + 12] + 0.8 * dp[m - 3]
            
            # The optimal strategy is to choose the one with the max probability.
            # Store it only if the probability is non-zero.
            res = max(prob_alpha, prob_beta)
            if res > 1e-9:  # Use a threshold for floating point
                new_dp[m] = res
        
        dp = new_dp

    # The final answer is the probability of success starting with M capital and n trades.
    final_probability = dp[M]

    # --- Output the final equation as requested ---
    print("--- Optimal Strategy Calculation ---")
    print(f"Starting with M=£{M} and n={n} trades. Goal is to reach £{target_money}.")
    
    # Values from the dp table for n-1 trades remaining
    val_alpha_succ = dp_n_minus_1[M + 1]
    val_alpha_fail = dp_n_minus_1[M - 1]
    val_beta_succ = dp_n_minus_1[M + 12]
    val_beta_fail = dp_n_minus_1[M - 3]
    
    prob_alpha_final = 0.0
    if M >= 1:
        prob_alpha_final = 0.6 * val_alpha_succ + 0.4 * val_alpha_fail
        print(f"\nStrategy Alpha (cost £1):")
        print(f"  - Prob. of success after gain (at £{M+1}): {val_alpha_succ:.6f}")
        print(f"  - Prob. of success after loss (at £{M-1}): {val_alpha_fail:.6f}")
        print(f"  - Resulting P(Alpha) = 0.60 * {val_alpha_succ:.6f} + 0.40 * {val_alpha_fail:.6f} = {prob_alpha_final:.6f}")
    else:
        print("\nStrategy Alpha is not affordable.")

    prob_beta_final = 0.0
    if M >= 3:
        prob_beta_final = 0.2 * val_beta_succ + 0.8 * val_beta_fail
        print(f"\nStrategy Beta (cost £3):")
        print(f"  - Prob. of success after gain (at £{M+12}): {val_beta_succ:.6f}")
        print(f"  - Prob. of success after loss (at £{M-3}):  {val_beta_fail:.6f}")
        print(f"  - Resulting P(Beta)  = 0.20 * {val_beta_succ:.6f} + 0.80 * {val_beta_fail:.6f} = {prob_beta_final:.6f}")
    else:
        print("\nStrategy Beta is not affordable.")
    
    print("\n--- Final Result ---")
    print(f"Optimal decision is to choose the strategy with the maximum probability.")
    print(f"Max Probability = max({prob_alpha_final:.6f}, {prob_beta_final:.6f})")
    print(f"The maximum probability of reaching £{target_money} is: {final_probability:.6f}")


# --- User-configurable parameters ---
# Initial Investment
INITIAL_INVESTMENT_M = 50
# Number of Trades
NUM_TRADES_N = 10
# --- End of parameters ---

solve_trading_problem(INITIAL_INVESTMENT_M, NUM_TRADES_N)
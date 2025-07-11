import sys

# It's possible for the recursion limit to be hit with very large n, though unlikely
# with this DP approach. Setting a higher limit is a safeguard.
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

    # Define probabilities and net changes for each strategy outcome
    prob_alpha_success = 0.60
    change_alpha_success = 1  # Cost £1, return £2 -> Net +£1
    prob_alpha_fail = 0.40
    change_alpha_fail = -1 # Cost £1, return £0 -> Net -£1

    prob_beta_success = 0.20
    change_beta_success = 12 # Cost £3, return £15 -> Net +£12
    prob_beta_fail = 0.80
    change_beta_fail = -3  # Cost £3, return £0 -> Net -£3

    # Determine the required size of the DP table for money (m)
    target_m = 2 * M
    # The max possible money is the starting amount plus the max possible gain over n trades.
    # We also need to be able to index the target amount, 2M.
    max_m_from_gains = M + 12 * n
    # Add a small buffer to prevent index out of bounds on m+12 access
    max_m = max(target_m, max_m_from_gains) + 13

    # dp[l][m]: max probability of success with m money after l trades.
    dp = [[0.0 for _ in range(max_m + 1)] for _ in range(n + 1)]

    # Base case: After n trades (l=n), success is 1.0 only if m is exactly 2M.
    if target_m <= max_m:
        dp[n][target_m] = 1.0

    # Fill the DP table by iterating backwards from the second to last trade
    for l in range(n - 1, -1, -1):
        for m in range(max_m):  # Iterate up to max_m-1, as we access m+1 etc.
            # --- Probability if choosing Strategy Alpha ---
            prob_alpha = 0.0
            if m >= 1:  # Must have £1 to invest
                # Get probability from the next state after a successful Alpha trade
                p_alpha_succ = dp[l + 1][m + change_alpha_success]
                # Get probability from the next state after a failed Alpha trade
                p_alpha_fail = dp[l + 1][m + change_alpha_fail]
                prob_alpha = prob_alpha_success * p_alpha_succ + prob_alpha_fail * p_alpha_fail

            # --- Probability if choosing Strategy Beta ---
            prob_beta = 0.0
            if m >= 3:  # Must have £3 to invest
                # Get probability from the next state after a successful Beta trade
                p_beta_succ = dp[l + 1][m + change_beta_success]
                # Get probability from the next state after a failed Beta trade
                p_beta_fail = dp[l + 1][m + change_beta_fail]
                prob_beta = prob_beta_success * p_beta_succ + prob_beta_fail * p_beta_fail

            # The optimal strategy is the one with the maximum probability
            dp[l][m] = max(prob_alpha, prob_beta)

    # The final answer is the probability at the start (l=0, m=M)
    final_prob = dp[0][M]

    # --- Output the results and the equation for the first move ---
    print(f"Solving for M = £{M} and n = {n} trades.")
    print(f"Goal: Reach exactly £{target_m} at the end.")
    print("-" * 50)
    print("Optimal Probability Calculation for the First Trade (l=0, m={M}):")
    print("T[m, l] = Max probability of success with £m after l trades.\n")

    # Values for Strategy Alpha at the start
    p_alpha_succ_start = dp[1][M + change_alpha_success] if M>=1 else 0
    p_alpha_fail_start = dp[1][M + change_alpha_fail] if M>=1 else 0
    prob_alpha_start = prob_alpha_success * p_alpha_succ_start + prob_alpha_fail * p_alpha_fail_start if M>=1 else 0

    print("If choosing Strategy Alpha:")
    print(f"  Prob = {prob_alpha_success:.2f} * T(m={M + change_alpha_success}, l=1) + {prob_alpha_fail:.2f} * T(m={M + change_alpha_fail}, l=1)")
    print(f"       = {prob_alpha_success:.2f} * {p_alpha_succ_start:.6f} + {prob_alpha_fail:.2f} * {p_alpha_fail_start:.6f}")
    print(f"       = {prob_alpha_start:.6f}")

    # Values for Strategy Beta at the start
    p_beta_succ_start = dp[1][M + change_beta_success] if M>=3 else 0
    p_beta_fail_start = dp[1][M + change_beta_fail] if M>=3 else 0
    prob_beta_start = prob_beta_success * p_beta_succ_start + prob_beta_fail * p_beta_fail_start if M>=3 else 0
    
    print("\nIf choosing Strategy Beta:")
    print(f"  Prob = {prob_beta_success:.2f} * T(m={M + change_beta_success}, l=1) + {prob_beta_fail:.2f} * T(m={M + change_beta_fail}, l=1)")
    print(f"       = {prob_beta_success:.2f} * {p_beta_succ_start:.6f} + {prob_beta_fail:.2f} * {p_beta_fail_start:.6f}")
    print(f"       = {prob_beta_start:.6f}")
    
    print("-" * 50)
    
    # Announce the optimal choice
    if prob_alpha_start > prob_beta_start:
        print("Optimal first move is Strategy Alpha.")
    elif prob_beta_start > prob_alpha_start:
        print("Optimal first move is Strategy Beta.")
    else:
        if final_prob > 0:
            print("Both strategies yield the same optimal probability. Either can be chosen.")
        else:
             print("It's impossible to guarantee reaching the target. Optimal probability is 0.")


    print(f"\nThe maximum probability of doubling the investment is: {final_prob:.6f}")


# --- Example Execution ---
# You can change these values to test different scenarios
initial_investment_M = 25
number_of_trades_n = 10
solve_trading_probability(initial_investment_M, number_of_trades_n)
import sys

def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    by executing exactly n trades.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        None. Prints the result to the console.
    """
    # The problem requires exactly n trades, and the goal is to reach exactly 2M.
    # We will use dynamic programming with a backward recursion approach.
    
    # Let dp[l][m] be the maximum probability of reaching the target (2M),
    # given that we have 'l' trades remaining and 'm' GBP.

    # The maximum possible money we can have: start with M, and n trades with max profit of 12 each.
    MAX_MONEY = M + 12 * n
    TARGET_MONEY = 2 * M

    # Check if the target is even reachable in principle.
    # The net change from any trade is an integer. So if M and 2M have different parities
    # after n trades, it's impossible. Net change per trade is M_final - M_initial.
    # Alpha: +1 or -1 (odd change). Beta: +12 or -3 (odd change).
    # After n trades, total change is sum of n odd numbers.
    # If n is even, sum is even. If n is odd, sum is odd.
    # Change needed = 2M - M = M.
    # If (n is even and M is odd) or (n is odd and M is even), it's impossible.
    if (n % 2) != (M % 2):
        print(f"It's impossible to reach the target. The parity of the required gain ({M}) does not match the parity of the number of trades ({n}).")
        print("The maximum probability is: 0.0")
        return
        
    # The DP table size: (n+1) rows for trades, (MAX_MONEY+1) columns for money.
    # Initialized to 0.0
    try:
        dp = [[0.0 for _ in range(MAX_MONEY + 1)] for _ in range(n + 1)]
    except MemoryError:
        print(f"Error: The required memory for the DP table is too large (n={n}, M={M}).", file=sys.stderr)
        print("Please try with smaller values.", file=sys.stderr)
        return

    # Base Case: l = 0 trades remaining.
    # Success is 1.0 only if we have exactly TARGET_MONEY.
    if 0 <= TARGET_MONEY <= MAX_MONEY:
        dp[0][TARGET_MONEY] = 1.0

    # Fill the DP table bottom-up (from l=1 to l=n).
    for l in range(1, n + 1):  # l = number of trades remaining
        for m in range(MAX_MONEY + 1):  # m = current money
            
            # --- Option 1: Use Strategy Alpha ---
            prob_alpha = 0.0
            if m >= 1:  # Must have at least £1 for the fee
                # Future states after one Alpha trade
                m_success = m + 1
                m_failure = m - 1
                prob_alpha = 0.6 * dp[l-1][m_success] + 0.4 * dp[l-1][m_failure]

            # --- Option 2: Use Strategy Beta ---
            prob_beta = 0.0
            if m >= 3:  # Must have at least £3 for the fee
                # Future states after one Beta trade
                m_success = m + 12
                m_failure = m - 3
                # Ensure the success state is within our table bounds
                prob_from_success = dp[l-1][m_success] if m_success <= MAX_MONEY else 0.0
                prob_beta = 0.2 * prob_from_success + 0.8 * dp[l-1][m_failure]
            
            # The optimal strategy is to choose the one with the maximum probability of success.
            dp[l][m] = max(prob_alpha, prob_beta)

    # The final answer is the state with n trades remaining and M initial capital.
    final_probability = dp[n][M]

    # --- Print the results and the final equation as requested ---
    print(f"Initial Investment M: £{M}")
    print(f"Number of Trades n: {n}")
    print(f"Target Investment 2M: £{TARGET_MONEY}")
    print("-" * 30)

    # Re-calculate the final step to show the equation
    final_prob_alpha = 0.0
    if M >= 1:
        final_prob_alpha = 0.6 * dp[n-1][M+1] + 0.4 * dp[n-1][M-1]

    final_prob_beta = 0.0
    if M >= 3:
        prob_from_success = dp[n-1][M+12] if M + 12 <= MAX_MONEY else 0.0
        final_prob_beta = 0.2 * prob_from_success + 0.8 * dp[n-1][M-3]

    print("Final Equation Breakdown (Optimal choice for the first trade):")
    print(f"T({n}, {M}) = max(Prob_Alpha, Prob_Beta)\n")

    print("1. Strategy Alpha:")
    if M >= 1:
        print(f"   Prob_Alpha = 0.60 * T({n-1}, {M+1}) + 0.40 * T({n-1}, {M-1})")
        print(f"   Prob_Alpha = 0.60 * {dp[n-1][M+1]:.5f} + 0.40 * {dp[n-1][M-1]:.5f}")
        print(f"   Resulting Probability = {final_prob_alpha:.5f}")
    else:
        print("   Cannot afford Strategy Alpha (Cost: £1)")

    print("\n2. Strategy Beta:")
    if M >= 3:
        print(f"   Prob_Beta = 0.20 * T({n-1}, {M+12}) + 0.80 * T({n-1}, {M-3})")
        prob_from_success_val = dp[n-1][M+12] if M+12 <= MAX_MONEY else 0.0
        print(f"   Prob_Beta = 0.20 * {prob_from_success_val:.5f} + 0.80 * {dp[n-1][M-3]:.5f}")
        print(f"   Resulting Probability = {final_prob_beta:.5f}")
    else:
        print("   Cannot afford Strategy Beta (Cost: £3)")

    print("-" * 30)
    print(f"The maximum probability of achieving the target is: {final_probability:.5f}")


# --- Example Usage ---
# You can change these values to test different scenarios.
INITIAL_INVESTMENT = 25
NUMBER_OF_TRADES = 10

solve_trading_probability(INITIAL_INVESTMENT, NUMBER_OF_TRADES)

# The time complexity of this algorithm is determined by the two nested loops.
# The outer loop runs 'n' times (for trades remaining from 1 to n).
# The inner loop runs 'MAX_MONEY + 1' times, where MAX_MONEY = M + 12*n.
# Total operations are proportional to n * (M + 12*n + 1).
# Therefore, the time complexity is O(n * (M + n)).
# Time Complexity: O(n*M + n^2)
<<<O(n(M+n))>>>
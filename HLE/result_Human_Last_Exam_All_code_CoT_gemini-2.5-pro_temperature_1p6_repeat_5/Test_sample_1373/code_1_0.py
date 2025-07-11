import sys

# It's recommended to increase recursion limit for larger N, although this DP is iterative.
# sys.setrecursionlimit(2000) 

def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        float: The maximum probability of reaching exactly 2*M GBP.
    """
    # The maximum possible gain from a single trade is +£12 (Strategy Beta).
    # This determines the upper bound for the money we can have.
    M_max = M + 12 * n
    
    # Target amount of money
    target_M = 2 * M

    # DP table: dp[l][m] is the max probability of success
    # given we are at trade l with m pounds.
    # Dimensions: (n+1) for trades 0..n, (M_max+1) for money 0..M_max.
    dp = [[0.0 for _ in range(M_max + 1)] for _ in range(n + 1)]

    # --- Base Case ---
    # At the end of n trades (l=n), success is having exactly target_M.
    # The target must be a reachable amount of money.
    if target_M <= M_max:
        dp[n][target_M] = 1.0

    # --- Iterative Step ---
    # Fill the DP table by iterating backwards from trade n-1 to 0.
    for l in range(n - 1, -1, -1):
        # Iterate through all possible money amounts at this stage.
        # An optimization could be to only iterate over reachable money states,
        # but this is simpler and correct.
        for m in range(M_max + 1):
            
            prob_alpha = 0.0
            prob_beta = 0.0

            # --- Strategy Alpha ---
            # Requires £1 fee. Money changes: m -> m+1 (60%) or m -> m-1 (40%)
            if m >= 1:
                # Check for array bounds, although our M_max setup should prevent this.
                if (m + 1 <= M_max) and (m - 1 >= 0):
                    prob_success_alpha = 0.6 * dp[l + 1][m + 1]
                    prob_failure_alpha = 0.4 * dp[l + 1][m - 1]
                    prob_alpha = prob_success_alpha + prob_failure_alpha

            # --- Strategy Beta ---
            # Requires £3 fee. Money changes: m -> m+12 (20%) or m -> m-3 (80%)
            if m >= 3:
                # Check for array bounds.
                if (m + 12 <= M_max) and (m - 3 >= 0):
                    prob_success_beta = 0.2 * dp[l + 1][m + 12]
                    prob_failure_beta = 0.8 * dp[l + 1][m - 3]
                    prob_beta = prob_success_beta + prob_failure_beta
            
            # The optimal choice is the one with the maximum probability of success.
            dp[l][m] = max(prob_alpha, prob_beta)

    # The final answer is the probability at the start (trade 0, money M).
    final_probability = dp[0][M]
    
    print(f"Initial investment M: {M} GBP")
    print(f"Number of trades n: {n}")
    print(f"Target investment 2M: {target_M} GBP")
    print(f"The optimal probability of success is: {final_probability}")


if __name__ == '__main__':
    # --- Parameters ---
    # The initial amount of money in GBP.
    initial_investment_M = 25
    # The exact number of trades to perform.
    number_of_trades_n = 10
    
    solve_trading_probability(initial_investment_M, number_of_trades_n)
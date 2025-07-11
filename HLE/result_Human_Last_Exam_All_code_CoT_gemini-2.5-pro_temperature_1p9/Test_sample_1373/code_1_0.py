import sys

def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    to 2M in exactly n trades using dynamic programming.
    
    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.
    """
    if M <= 0 or n < 0:
        print("Initial investment and number of trades must be non-negative.")
        return

    target = 2 * M
    # The maximum possible amount of money occurs if we always pick Beta and win.
    max_money_from_trades = M + 12 * n
    
    # The DP table for money needs to accommodate the target and any potential future states.
    # A state (l, m) could look up m+12 in the next step. So the bound needs to accommodate that.
    money_bound = max(target, max_money_from_trades)
    
    # Initialize DP table: dp[l][m]
    # l: trades completed (0 to n)
    # m: money held (0 to money_bound + 12, to be safe for m+12 lookups)
    dp_table_cols = money_bound + 13
    dp = [[0.0 for _ in range(dp_table_cols)] for _ in range(n + 1)]

    # Base Case: At the end of n trades (l=n), the probability of success is 1.0
    # if we have exactly 2M, and 0 otherwise.
    if target < dp_table_cols:
        dp[n][target] = 1.0

    # Iterate backwards from the second to last trade to the first trade
    for l in range(n - 1, -1, -1):
        # Iterate over all possible money amounts at this stage
        # The max money possible after l trades is M + 12*l.
        # We loop a bit higher to ensure all dependencies are met.
        current_loop_bound = M + 12 * l + 1 
        for m in range(current_loop_bound):
            # Calculate probability of success if we choose Strategy Alpha
            p_alpha = 0.0
            if m >= 1:
                prob_success_alpha = dp[l + 1][m + 1]
                prob_failure_alpha = dp[l + 1][m - 1]
                p_alpha = 0.6 * prob_success_alpha + 0.4 * prob_failure_alpha
            
            # Calculate probability of success if we choose Strategy Beta
            p_beta = 0.0
            if m >= 3:
                prob_success_beta = dp[l + 1][m + 12]
                prob_failure_beta = dp[l + 1][m - 3]
                p_beta = 0.2 * prob_success_beta + 0.8 * prob_failure_beta
            
            # Choose the optimal strategy
            if m >= 3:
                dp[l][m] = max(p_alpha, p_beta)
            elif m >= 1:
                dp[l][m] = p_alpha
            else: # m < 1, cannot trade
                dp[l][m] = 0.0
                
    # The final answer is the probability at the initial state (0 trades, M money)
    final_probability = dp[0][M]
    
    print(f"Initial Investment M = {M} GBP")
    print(f"Number of Trades n = {n}")
    print(f"Target Investment = {target} GBP")
    print("-" * 30)

    print("Analysis for the first trade decision:")
    
    # Recalculate and display the choice for the first trade for clarity
    p_alpha_initial = 0.0
    if M >= 1:
        prob_success_alpha = dp[1][M + 1]
        prob_failure_alpha = dp[1][M - 1]
        p_alpha_initial = 0.6 * prob_success_alpha + 0.4 * prob_failure_alpha
        print(f"If choosing Alpha (cost £1):")
        print(f"  Success changes money to £{M+1}, failure to £{M-1}.")
        print(f"  P(success from start with Alpha) = 0.60 * P(success | £{M+1}, 1 trade) + 0.40 * P(success | £{M-1}, 1 trade)")
        print(f"  P(success with Alpha) = 0.60 * {prob_success_alpha:.4f} + 0.40 * {prob_failure_alpha:.4f} = {p_alpha_initial:.4f}")

    p_beta_initial = 0.0
    if M >= 3:
        prob_success_beta = dp[1][M + 12]
        prob_failure_beta = dp[1][M - 3]
        p_beta_initial = 0.2 * prob_success_beta + 0.8 * prob_failure_beta
        print(f"If choosing Beta (cost £3):")
        print(f"  Success changes money to £{M+12}, failure to £{M-3}.")
        print(f"  P(success from start with Beta) = 0.20 * P(success | £{M+12}, 1 trade) + 0.80 * P(success | £{M-3}, 1 trade)")
        print(f"  P(success with Beta) = 0.20 * {prob_success_beta:.4f} + 0.80 * {prob_failure_beta:.4f} = {p_beta_initial:.4f}")

    print("-" * 30)
    print(f"The maximum probability of achieving exactly £{target} is: {final_probability:.4f}")


if __name__ == '__main__':
    # You can change these values to test different scenarios
    initial_money = 25
    num_trades = 4
    solve_trading_probability(initial_money, num_trades)

    # A scenario that is likely impossible (should result in 0)
    # print("\n" + "="*50 + "\n")
    # solve_trading_probability(M=100, n=2)

    # A scenario that might have a solution
    # print("\n" + "="*50 + "\n")
    # solve_trading_probability(M=19, n=1) # Should be 0.6

# The complexity is determined by the size of the DP table and the work per cell.
# The table size is (n+1) rows and O(max(M, M+n)) = O(M+n) columns.
# Work per cell is O(1).
# Total time complexity = O(n * (M+n)) = O(n*M + n^2)
<<<O(n*M + n^2)>>>
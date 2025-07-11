import collections

def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        float: The maximum probability of reaching exactly 2*M GBP.
    """
    # The maximum profit from one trade is from Beta: 15 (return) - 3 (fee) = 12
    # A safe upper bound for money is the initial amount plus max profit from n trades.
    max_money = M + n * 12
    target_money = 2 * M

    # dp[i][j]: max probability of success with i trades left and j pounds.
    # We use two rows to optimize space: dp_prev for i-1, dp_curr for i.
    dp_prev = collections.defaultdict(float)
    
    # Base case: 0 trades left.
    # Success only if money is exactly the target.
    if target_money <= max_money:
        dp_prev[target_money] = 1.0

    # Iterate from 1 to n trades remaining
    for i in range(1, n + 1):
        dp_curr = collections.defaultdict(float)
        # Iterate through all possible money amounts
        for m in range(max_money + 1):
            # Calculate probability of success if we choose Strategy Alpha
            prob_alpha = 0.0
            if m >= 1:
                # 60% chance of success (money becomes m+1)
                # 40% chance of failure (money becomes m-1)
                prob_alpha = 0.6 * dp_prev[m + 1] + 0.4 * dp_prev[m - 1]

            # Calculate probability of success if we choose Strategy Beta
            prob_beta = 0.0
            if m >= 3:
                # 20% chance of success (money becomes m+12)
                # 80% chance of failure (money becomes m-3)
                prob_beta = 0.2 * dp_prev[m + 12] + 0.8 * dp_prev[m - 3]
            
            # The optimal strategy is to choose the one with the max probability
            dp_curr[m] = max(prob_alpha, prob_beta)
        
        dp_prev = dp_curr

    # The final answer is the probability of success starting with M pounds and n trades
    final_probability = dp_prev[M]
    
    print(f"Initial Investment (M): £{M}")
    print(f"Number of Trades (n): {n}")
    print(f"Target Investment (2M): £{target_money}")
    print(f"The maximum probability of reaching the target is: {final_probability}")

if __name__ == '__main__':
    # Example values for M and n
    initial_investment = 25
    num_trades = 5
    solve_trading_problem(initial_investment, num_trades)

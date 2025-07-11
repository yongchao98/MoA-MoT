import sys

# It's possible to run into recursion limits for very large n, though the iterative DP approach avoids this.
# For a recursive solution (memoization), one might need to increase this.
# sys.setrecursionlimit(2000) 

def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades.

    Args:
        M (int): The initial investment in GBP. Must be an integer > 0.
        n (int): The exact number of trades to be executed. Must be an integer >= 0.

    Returns:
        float: The maximum probability of success.
    """
    if M <= 0 or n < 0:
        print("Initial investment M and number of trades n must be non-negative.")
        return 0.0

    if n == 0:
        return 1.0 if M == 2 * M else 0.0 # 1.0 only if M is 0, which is disallowed. Effectively 0.

    # Constants for the trading strategies
    ALPHA_FEE = 1
    ALPHA_SUCCESS_PROB = 0.6
    ALPHA_NET_PROFIT_SUCCESS = 1  # Net is 2(return) - 1(fee) = 1
    ALPHA_NET_PROFIT_FAIL = -1 # Net is 0(return) - 1(fee) = -1

    BETA_FEE = 3
    BETA_SUCCESS_PROB = 0.2
    BETA_NET_PROFIT_SUCCESS = 12 # Net is 15(return) - 3(fee) = 12
    BETA_NET_PROFIT_FAIL = -3  # Net is 0(return) - 3(fee) = -3
    
    # The target amount of money
    TARGET_MONEY = 2 * M

    # Determine the maximum money we might need to track in our DP table.
    # At any step i (trades left), the min money we can have to reach target is 2M - 12*i
    # The max money is 2M + 3*i. So we need to consider states up to 2M + 3*n.
    # Any money amount above this is not on an optimal path back to 2M.
    max_money = TARGET_MONEY + 3 * n
    
    # dp[i][j]: max probability of success with i trades left and j money.
    # Dimensions: (n+1) trades, (max_money+1) money
    dp = [[0.0 for _ in range(max_money + 1)] for _ in range(n + 1)]

    # Base case: i = 0 trades left.
    # Probability is 1.0 if we are exactly at the target, 0.0 otherwise.
    if TARGET_MONEY <= max_money:
        dp[0][TARGET_MONEY] = 1.0

    # Fill the DP table iteratively
    # i is the number of trades remaining
    for i in range(1, n + 1):
        # j is the current amount of money
        for j in range(max_money + 1):
            prob_alpha = 0.0
            prob_beta = 0.0

            # Option 1: Strategy Alpha
            if j >= ALPHA_FEE:
                # Money after success
                next_money_s = j + ALPHA_NET_PROFIT_SUCCESS
                # Money after failure
                next_money_f = j + ALPHA_NET_PROFIT_FAIL

                # Get probabilities from the previous state (i-1 trades left)
                # Check boundaries to avoid index errors
                p_s_alpha = dp[i-1][next_money_s] if 0 <= next_money_s <= max_money else 0.0
                p_f_alpha = dp[i-1][next_money_f] if 0 <= next_money_f <= max_money else 0.0

                prob_alpha = ALPHA_SUCCESS_PROB * p_s_alpha + (1 - ALPHA_SUCCESS_PROB) * p_f_alpha

            # Option 2: Strategy Beta
            if j >= BETA_FEE:
                # Money after success
                next_money_s = j + BETA_NET_PROFIT_SUCCESS
                # Money after failure
                next_money_f = j + BETA_NET_PROFIT_FAIL
                
                p_s_beta = dp[i-1][next_money_s] if 0 <= next_money_s <= max_money else 0.0
                p_f_beta = dp[i-1][next_money_f] if 0 <= next_money_f <= max_money else 0.0

                prob_beta = BETA_SUCCESS_PROB * p_s_beta + (1 - BETA_SUCCESS_PROB) * p_f_beta

            # The optimal strategy maximizes the probability
            dp[i][j] = max(prob_alpha, prob_beta)

    # The final answer is the probability starting with M money and n trades
    result_prob = dp[n][M]
    return result_prob

if __name__ == '__main__':
    # Example values for M and n
    # For M=3, n=3, the optimal strategy involves always choosing Alpha.
    # The only way to get a net profit of 3 is 3 successful Alpha trades (+1+1+1=3).
    # Path: 3 ->(A)-> 4 ->(A)-> 5 ->(A)-> 6.
    # The probability of this specific path is 0.6 * 0.6 * 0.6 = 0.216.
    # DP will confirm if this is the optimal path.
    initial_investment_M = 3
    num_trades_n = 3

    probability = solve_trading_probability(initial_investment_M, num_trades_n)

    print(f"Starting with £{initial_investment_M} and executing {num_trades_n} trades.")
    print(f"The target is to reach exactly £{2 * initial_investment_M}.")
    print(f"The maximum probability of achieving this target is: {probability:.10f}")

    # Another example where the probability should be non-zero
    # M=1, n=1 -> Net profit of +1 is needed.
    # One successful Alpha trade gives +1. Fee is 1, affordable. Prob = 0.6
    print("\n--- Another Example ---")
    initial_investment_M_2 = 1
    num_trades_n_2 = 1
    probability_2 = solve_trading_probability(initial_investment_M_2, num_trades_n_2)
    print(f"Starting with £{initial_investment_M_2} and executing {num_trades_n_2} trades.")
    print(f"The target is to reach exactly £{2 * initial_investment_M_2}.")
    print(f"The maximum probability of achieving this target is: {probability_2:.10f}")


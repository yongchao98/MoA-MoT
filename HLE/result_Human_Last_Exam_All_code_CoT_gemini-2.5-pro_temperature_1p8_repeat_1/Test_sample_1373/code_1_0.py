import sys

# It is recommended to increase the recursion limit for larger n,
# although the iterative DP approach avoids this need.
# sys.setrecursionlimit(2000)

def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        float: The maximum probability of success (reaching 2*M).
    """
    target_money = 2 * M
    # The maximum possible money occurs if we start with M and always have
    # successful Beta trades (net profit of 12 each).
    max_money = M + 12 * n

    # dp_table[l][m] stores the max probability of success with l trades
    # remaining and m GBP.
    dp_table = [[0.0 for _ in range(max_money + 1)] for _ in range(n + 1)]

    # Base case: l = 0 trades remaining.
    # We succeed only if our money is exactly the target amount.
    if 0 <= target_money <= max_money:
        dp_table[0][target_money] = 1.0

    # Fill the DP table using backward induction from l=1 to n.
    for l in range(1, n + 1):  # l = number of trades remaining
        for m in range(max_money + 1):  # m = current money
            
            # --- Calculate probability of success for Strategy Alpha ---
            prob_alpha = 0.0
            if m >= 1: # Check if we can afford the fee for Alpha
                # T[l-1][m+1] and T[l-1][m-1] are the probabilities of success
                # from the states we would land in after the Alpha trade.
                p_success = dp_table[l-1][m + 1] if (m + 1) <= max_money else 0.0
                p_fail = dp_table[l-1][m - 1] if (m - 1) >= 0 else 0.0
                prob_alpha = 0.6 * p_success + 0.4 * p_fail

            # --- Calculate probability of success for Strategy Beta ---
            prob_beta = 0.0
            if m >= 3: # Check if we can afford the fee for Beta
                p_success = dp_table[l-1][m + 12] if (m + 12) <= max_money else 0.0
                p_fail = dp_table[l-1][m - 3] if (m - 3) >= 0 else 0.0
                prob_beta = 0.2 * p_success + 0.8 * p_fail

            # --- Determine the optimal strategy and store its probability ---
            # We must choose the strategy that maximizes our chance of success.
            if m >= 3:
                dp_table[l][m] = max(prob_alpha, prob_beta)
            elif m >= 1:
                dp_table[l][m] = prob_alpha
            else: # m < 1, cannot afford any trade, so probability of success is 0.
                dp_table[l][m] = 0.0

    # The final answer is the probability of success starting with M money and n trades.
    final_probability = dp_table[n][M]
    return final_probability

if __name__ == '__main__':
    # Define the initial investment and number of trades.
    # You can change these values to test other scenarios.
    initial_investment_M = 25
    num_trades_n = 10

    print(f"Calculating the max probability of doubling an initial investment of £{initial_investment_M}")
    print(f"by executing exactly {num_trades_n} trades.\n")
    
    # Calculate the probability
    max_prob = solve_trading_probability(initial_investment_M, num_trades_n)
    
    # Output the final result. As the problem asks, the numbers involved in
    # the final result's calculation are implicitly handled within the code logic.
    # The final equation is T[n][M] = result, where T is our DP table.
    # Here, T[10][25] is the result we print.
    print(f"Initial State: (Trades remaining={num_trades_n}, Money=£{initial_investment_M})")
    print(f"Target Money: £{2 * initial_investment_M}")
    print(f"Maximum Probability of Success: {max_prob:.6f}")

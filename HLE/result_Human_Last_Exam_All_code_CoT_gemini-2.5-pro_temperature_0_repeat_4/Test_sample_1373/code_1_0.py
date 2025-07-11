import numpy as np

def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    by executing exactly n trades.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        float: The maximum probability of achieving the target of 2M.
    """
    # The maximum possible amount of money after n trades
    max_money = M + 12 * n
    target_money = 2 * M

    # DP table: T[m][l] is the max probability of success with £m after l trades.
    # Dimensions: (max_money + 1) rows for money, (n + 1) columns for trades.
    T = np.zeros((max_money + 1, n + 1))

    # Base case: After n trades (l=n)
    # Success is 1.0 only if money is exactly 2M, otherwise 0.
    if 0 <= target_money <= max_money:
        T[target_money][n] = 1.0

    # Fill the DP table by iterating backwards from l = n-1 to 0
    for l in range(n - 1, -1, -1):
        for m in range(max_money + 1):
            # Probability of success if choosing Strategy Alpha
            p_alpha = 0.0
            if m >= 1:
                # Net profit is £1 (cost £1, return £2)
                # Net loss is £1 (cost £1, return £0)
                # Ensure indices are within bounds
                p_success_alpha = T[m + 1][l + 1] if (m + 1) <= max_money else 0.0
                p_fail_alpha = T[m - 1][l + 1] if (m - 1) >= 0 else 0.0
                p_alpha = 0.6 * p_success_alpha + 0.4 * p_fail_alpha

            # Probability of success if choosing Strategy Beta
            p_beta = 0.0
            if m >= 3:
                # Net profit is £12 (cost £3, return £15)
                # Net loss is £3 (cost £3, return £0)
                # Ensure indices are within bounds
                p_success_beta = T[m + 12][l + 1] if (m + 12) <= max_money else 0.0
                p_fail_beta = T[m - 3][l + 1] if (m - 3) >= 0 else 0.0
                p_beta = 0.2 * p_success_beta + 0.8 * p_fail_beta

            # The optimal strategy maximizes the probability of success
            T[m][l] = max(p_alpha, p_beta)

    # The final answer is the probability at the start (M money, 0 trades)
    return T[M][0]

if __name__ == '__main__':
    # --- Problem Parameters ---
    # Initial Investment in GBP
    initial_investment_M = 25
    # Number of trades
    num_trades_n = 10
    # Target investment
    target_investment = 2 * initial_investment_M

    # --- Calculation ---
    max_probability = solve_trading_probability(initial_investment_M, num_trades_n)

    # --- Outputting the final equation's parameters and result ---
    print("--- Trading Problem Setup ---")
    print(f"Initial Investment (M): £{initial_investment_M}")
    print(f"Number of Trades (n): {num_trades_n}")
    print(f"Target Investment (2M): £{target_investment}")
    print("\n--- Optimal Strategy Result ---")
    print(f"The maximum probability of reaching the target is: {max_probability:.6f}")

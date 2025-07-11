import sys

def calculate_max_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        float: The maximum probability of success.
    """
    if M <= 0 or n < 0:
        return 0.0

    # The maximum net gain per trade is +12 from a successful Beta trade.
    max_gain_per_trade = 12
    # The maximum possible money we might need to track in our DP table.
    # It must be large enough to hold the target (2*M) and any possible
    # capital state (M + n * max_gain).
    max_money = max(2 * M, M + n * max_gain_per_trade)

    # dp[i][j]: max probability of reaching 2M at the end of n trades,
    # given that 'i' trades are completed and we have 'j' GBP.
    # Initialize DP table with zeros.
    dp = [[0.0] * (max_money + 1) for _ in range(n + 1)]

    # Base Case: At the end of n trades (i=n).
    # Success is 100% if we have exactly 2M, 0% otherwise.
    target_money = 2 * M
    if target_money <= max_money:
        dp[n][target_money] = 1.0

    # Iterate backwards from trade n-1 down to 0.
    for i in range(n - 1, -1, -1):
        # Iterate over all possible money amounts 'j' at step 'i'.
        # The maximum possible capital at step i is M + i * max_gain_per_trade.
        # We iterate up to max_money for simpler logic.
        for j in range(max_money + 1):
            
            # --- Calculate probability of success for Strategy Alpha ---
            p_alpha = 0.0
            fee_alpha = 1
            if j >= fee_alpha:
                # Capital after winning: j - fee + profit = j - 1 + 2 = j + 1
                # Capital after losing: j - fee + profit = j - 1 + 0 = j - 1
                capital_win = j + 1
                capital_lose = j - 1
                
                # Probability of success from the resulting states
                prob_from_win = dp[i + 1][capital_win] if capital_win <= max_money else 0.0
                prob_from_lose = dp[i + 1][capital_lose] if capital_lose >= 0 else 0.0
                
                p_alpha = 0.6 * prob_from_win + 0.4 * prob_from_lose

            # --- Calculate probability of success for Strategy Beta ---
            p_beta = 0.0
            fee_beta = 3
            if j >= fee_beta:
                # Capital after winning: j - fee + profit = j - 3 + 15 = j + 12
                # Capital after losing: j - fee + profit = j - 3 + 0 = j - 3
                capital_win = j + 12
                capital_lose = j - 3

                # Probability of success from the resulting states
                prob_from_win = dp[i + 1][capital_win] if capital_win <= max_money else 0.0
                prob_from_lose = dp[i + 1][capital_lose] if capital_lose >= 0 else 0.0

                p_beta = 0.2 * prob_from_win + 0.8 * prob_from_lose

            # Choose the optimal strategy for state (i, j)
            dp[i][j] = max(p_alpha, p_beta)

    # The final answer is the probability starting with M money and 0 trades.
    return dp[0][M]

if __name__ == '__main__':
    # You can change these values to test different scenarios
    initial_investment_M = 25
    number_of_trades_n = 10

    # Increase recursion limit for potentially deep state spaces, though our iterative
    # approach avoids this issue. It's good practice for DP problems in general.
    sys.setrecursionlimit(2000)

    print(f"Calculating for M = {initial_investment_M}, n = {number_of_trades_n}")
    
    final_probability = calculate_max_probability(initial_investment_M, number_of_trades_n)
    
    print(f"The maximum probability of doubling the investment is: {final_probability:.6f}")
    
    # Time complexity explanation
    # The algorithm iterates 'n' times in the outer loop.
    # The inner loop iterates up to 'max_money' times.
    # max_money is O(M + n).
    # Therefore, the total time complexity is O(n * (M + n)).
    complexity = "O(n * (M + n))"
    print(f"The time complexity of this algorithm is: {complexity}")

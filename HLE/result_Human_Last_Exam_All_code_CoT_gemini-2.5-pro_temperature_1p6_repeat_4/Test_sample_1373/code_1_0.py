import sys

# Increase recursion limit for potentially deep financial paths, although our DP is iterative.
# This is more of a safeguard.
sys.setrecursionlimit(2000)

def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    after exactly n trades using dynamic programming.
    """

    # We need a capital buffer for intermediate calculations, especially for Strategy Beta's large profit.
    # The highest capital that could rationally lead to the 2M target is around 2M + 12.
    # A safe max capital for our DP table considers the max possible jump from any state we compute.
    # Max profit per trade is 12. So, an n-trade path could at most increase capital by n*12.
    # We set a generous upper bound for the capital to handle all reachable states.
    MAX_CAPITAL = 2 * M + n * 12 + 15

    # dp[i][j]: max probability of success with i trades remaining and capital j
    try:
        dp = [[0.0 for _ in range(MAX_CAPITAL + 1)] for _ in range(n + 1)]
    except MemoryError:
        print(f"Error: The required DP table is too large for memory (M={M}, n={n}).")
        print("Please try with smaller values of M and n.")
        return

    # Base case: 0 trades remaining
    # Success only if capital is exactly 2*M.
    if 2 * M <= MAX_CAPITAL:
        dp[0][2 * M] = 1.0

    # Fill the DP table using backward induction
    for i in range(1, n + 1):  # i = trades remaining
        for j in range(MAX_CAPITAL + 1):  # j = current capital

            # Probability of success if choosing Strategy Alpha
            p_alpha = 0.0
            if j >= 1:
                # Calculate probability ensuring we don't go out of bounds
                prob_alpha_success = 0.60 * dp[i - 1][j + 1] if (j + 1 <= MAX_CAPITAL) else 0.0
                prob_alpha_fail = 0.40 * dp[i - 1][j - 1] if (j - 1 >= 0) else 0.0
                p_alpha = prob_alpha_success + prob_alpha_fail

            # Probability of success if choosing Strategy Beta
            p_beta = 0.0
            if j >= 3:
                # Calculate probability ensuring we don't go out of bounds
                prob_beta_success = 0.20 * dp[i - 1][j + 12] if (j + 12 <= MAX_CAPITAL) else 0.0
                prob_beta_fail = 0.80 * dp[i - 1][j - 3] if (j - 3 >= 0) else 0.0
                p_beta = prob_beta_success + prob_beta_fail

            # Determine the optimal strategy for state (i, j)
            # We must make a trade. If not affordable, probability of success is 0.
            if j < 1:
                dp[i][j] = 0.0
            elif j < 3:
                dp[i][j] = p_alpha
            else:
                dp[i][j] = max(p_alpha, p_beta)

    # The result is the probability starting with M capital and n trades.
    final_probability = dp[n][M]

    # --- Output the details of the final calculation as requested ---
    print(f"Goal: Start with £{M}, execute {n} trades, and end with exactly £{2*M}.\n")
    print("This requires calculating the optimal strategy for each possible state.")
    print("The final probability is found at T[M, n] in our table.")
    print("-" * 40)
    print(f"Final Step Calculation: From {n-1} trades remaining to {n} trades remaining")
    print(f"Starting state: Capital = £{M}, Trades Remaining = {n}")
    print("-" * 40)

    # Re-calculate the choice at the final step for printing
    if n > 0:
        val_alpha_win_state = dp[n - 1][M + 1]
        val_alpha_lose_state = dp[n - 1][M - 1]
        p_alpha_final = 0.60 * val_alpha_win_state + 0.40 * val_alpha_lose_state

        val_beta_win_state = dp[n - 1][M + 12]
        val_beta_lose_state = dp[n - 1][M - 3]
        p_beta_final = 0.20 * val_beta_win_state + 0.80 * val_beta_lose_state

        if M < 1:
            print("Initial capital is less than £1. Cannot afford any trades.")
            print("Probability of success = 0")
        elif M < 3:
            print("Only Strategy Alpha is affordable:")
            print(f"Prob(Alpha) = 0.60 * T[£{M+1}, {n-1}] + 0.40 * T[£{M-1}, {n-1}]")
            print(f"Prob(Alpha) = 0.60 * {val_alpha_win_state:.6f} + 0.40 * {val_alpha_lose_state:.6f}")
            print(f"Resulting Probability = {p_alpha_final:.6f}")
        else:
            print("Strategy Alpha Calculation:")
            print(f"Prob(Alpha) = 0.60 * T[£{M+1}, {n-1}] + 0.40 * T[£{M-1}, {n-1}]")
            print(f"Prob(Alpha) = 0.60 * {val_alpha_win_state:.6f} + 0.40 * {val_alpha_lose_state:.6f} = {p_alpha_final:.6f}")
            print("\nStrategy Beta Calculation:")
            print(f"Prob(Beta)  = 0.20 * T[£{M+12}, {n-1}] + 0.80 * T[£{M-3}, {n-1}]")
            print(f"Prob(Beta)  = 0.20 * {val_beta_win_state:.6f} + 0.80 * {val_beta_lose_state:.6f} = {p_beta_final:.6f}")
            print(f"\nOptimal choice is the maximum: max({p_alpha_final:.6f}, {p_beta_final:.6f})")

    print("-" * 40)
    print(f"Final Maximum Probability: {final_probability:.6f}")


if __name__ == '__main__':
    # Example values for the trading problem
    initial_investment_M = 25
    number_of_trades_n = 10
    solve_trading_probability(initial_investment_M, number_of_trades_n)
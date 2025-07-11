import sys

def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.

    The time complexity of this algorithm is O(n * (M + n)).
    - The state is defined by (money, trades_done), leading to a 2D DP table.
    - The number of trades 'n' is one dimension.
    - The range of money 'm' is the other. The maximum money, M_max, can be approximated
      by max(2*M, M + 12*n), which is on the order of O(M + n).
    - We iterate through each state once, performing O(1) work.
    - Total complexity = O(n * M_max) = O(n * (M + n)) = O(Mn + n^2).
    """
    # Handle the trivial case where no trades are to be made.
    if n == 0:
        # Success only if the initial amount is already the target (2*M), which implies M=0.
        if M == 0:
            print("With 0 trades and 0 initial GBP, the target of 0 is met. Equation: 1.0")
            return 1.0
        else:
            print(f"With 0 trades, it's impossible to turn {M} GBP into {2*M} GBP. Probability = 0.0")
            return 0.0

    M_target = 2 * M
    # Determine the maximum possible money amount to size the DP table.
    M_max = max(M_target, M + 12 * n)

    # Initialize the DP table T[m][l].
    # T[m][l] stores the max probability of success starting with money 'm' after 'l' trades.
    T = [[0.0 for _ in range(n + 1)] for _ in range(M_max + 1)]

    # Base Case: At the end of the n-th trade (l=n).
    if M_target <= M_max:
        T[M_target][n] = 1.0

    # DP Calculation: Iterate backwards from l=n-1 to l=0.
    for l in range(n - 1, -1, -1):
        for m in range(M_max + 1):
            # Probability of success for Strategy Alpha
            prob_alpha = 0.0
            if m >= 1:
                prob_win_alpha = T[m + 1][l + 1]
                prob_lose_alpha = T[m - 1][l + 1]
                prob_alpha = 0.6 * prob_win_alpha + 0.4 * prob_lose_alpha

            # Probability of success for Strategy Beta
            prob_beta = 0.0
            if m >= 3:
                # Check if the win outcome is within the table bounds
                prob_win_beta = T[m + 12][l + 1] if m + 12 <= M_max else 0.0
                prob_lose_beta = T[m - 3][l + 1]
                prob_beta = 0.2 * prob_win_beta + 0.8 * prob_lose_beta

            # Optimal strategy chooses the action with the maximum probability of success.
            T[m][l] = max(prob_alpha, prob_beta)

    final_probability = T[M][0]

    # --- Output ---
    print(f"Initial investment: M = {M} GBP, Number of trades: n = {n}")
    print(f"Target investment: 2M = {M_target} GBP")
    print(f"The maximum probability of achieving the target is: {final_probability:.6f}")
    
    # Explain the first optimal move and show the calculation as an equation.
    if final_probability > 0:
        prob_alpha_initial = 0.0
        if M >= 1:
            prob_alpha_initial = 0.6 * T[M + 1][1] + 0.4 * T[M - 1][1]
        
        prob_beta_initial = 0.0
        if M >= 3:
            prob_beta_initial = 0.2 * T[M + 12][1] + 0.8 * T[M - 3][1] if M + 12 <= M_max else 0.0

        print("\n--- Analysis of the First Trade ---")
        if prob_alpha_initial >= prob_beta_initial:
            print("Optimal first move: Strategy Alpha")
            print(f"Equation: 0.6 * P(success | {M+1} GBP, 1 trade) + 0.4 * P(success | {M-1} GBP, 1 trade)")
            print(f"Result: 0.6 * {T[M+1][1]:.4f} + 0.4 * {T[M-1][1]:.4f} = {prob_alpha_initial:.6f}")
        else:
            print("Optimal first move: Strategy Beta")
            print(f"Equation: 0.2 * P(success | {M+12} GBP, 1 trade) + 0.8 * P(success | {M-3} GBP, 1 trade)")
            print(f"Result: 0.2 * {T[M+12][1]:.4f} + 0.8 * {T[M-3][1]:.4f} = {prob_beta_initial:.6f}")

if __name__ == '__main__':
    # Example scenario: Start with M=25 GBP and execute n=5 trades.
    # The goal is to reach exactly 50 GBP.
    initial_investment_M = 25
    number_of_trades_n = 5
    solve_trading_problem(initial_investment_M, number_of_trades_n)
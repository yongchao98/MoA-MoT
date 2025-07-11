import sys

def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.
    """
    # Set a high recursion limit for potentially deep state dependencies, though our DP is iterative.
    # sys.setrecursionlimit(2000) 

    target_m = 2 * M
    
    # The maximum possible amount of money is the initial capital plus the max profit from n trades.
    # Max profit per trade is +12 from Strategy Beta.
    m_max = M + 12 * n

    # If the target is unreachable even with the best-case scenario, probability is 0.
    if target_m > m_max:
        print(f"It's impossible to reach {target_m} GBP from {M} GBP in {n} trades.")
        print("Maximum possible amount is {m_max} GBP.")
        print("The probability of success is: 0.0")
        return

    # DP table: T[l][m] stores the max probability of success
    # given 'l' trades are completed and we have 'm' money.
    # Dimensions: (n + 1) trades, (m_max + 1) money amounts.
    T = [[0.0 for _ in range(m_max + 1)] for _ in range(n + 1)]

    # Base case: After n trades (l=n), success is 1.0 only if money is exactly 2M.
    if target_m <= m_max:
        T[n][target_m] = 1.0

    # Fill the DP table by iterating backwards from l = n-1 to 0.
    for l in range(n - 1, -1, -1):
        # Iterate over all possible money amounts m for the current trade number l.
        for m in range(m_max + 1):
            
            # --- Strategy Alpha ---
            p_alpha = 0.0
            # Check if we can afford the £1 fee.
            if m >= 1:
                # Success: money becomes m - 1 (fee) + 2 (return) = m + 1
                # Failure: money becomes m - 1 (fee) + 0 (return) = m - 1
                prob_from_success = 0.0
                if m + 1 <= m_max:
                    prob_from_success = 0.60 * T[l + 1][m + 1]

                prob_from_failure = 0.0
                if m - 1 >= 0:
                    prob_from_failure = 0.40 * T[l + 1][m - 1]
                
                p_alpha = prob_from_success + prob_from_failure

            # --- Strategy Beta ---
            p_beta = 0.0
            # Check if we can afford the £3 fee.
            if m >= 3:
                # Success: money becomes m - 3 (fee) + 15 (return) = m + 12
                # Failure: money becomes m - 3 (fee) + 0 (return) = m - 3
                prob_from_success = 0.0
                if m + 12 <= m_max:
                    prob_from_success = 0.20 * T[l + 1][m + 12]

                prob_from_failure = 0.0
                if m - 3 >= 0:
                    prob_from_failure = 0.80 * T[l + 1][m - 3]
                
                p_beta = prob_from_success + prob_from_failure

            # The optimal strategy is the one that maximizes the probability.
            T[l][m] = max(p_alpha, p_beta)
            
    # The final answer is the probability at the start (0 trades done, M money).
    final_probability = T[0][M]
    
    print(f"Goal: Reach {target_m} GBP from {M} GBP in {n} trades.\n")
    
    # As requested, printing the numbers for the final equation (at the first step)
    print(f"--- Analysis of the First Trade (starting with {M} GBP) ---")
    
    # Recalculate probabilities for the initial state for detailed output
    p_alpha_initial = 0.0
    if M >= 1:
        p_alpha_initial = 0.60 * T[1][M + 1] + 0.40 * T[1][M - 1]
        print(f"Option Alpha: 0.60 * P(success | {M+1} GBP, {n-1} trades left) + 0.40 * P(success | {M-1} GBP, {n-1} trades left)")
        print(f"            = 0.60 * {T[1][M + 1]:.4f} + 0.40 * {T[1][M - 1]:.4f} = {p_alpha_initial:.4f}")
    else:
        print("Option Alpha: Not affordable.")

    p_beta_initial = 0.0
    if M >= 3:
        p_beta_initial = 0.20 * T[1][M + 12] + 0.80 * T[1][M - 3]
        print(f"Option Beta:  0.20 * P(success | {M+12} GBP, {n-1} trades left) + 0.80 * P(success | {M-3} GBP, {n-1} trades left)")
        print(f"            = 0.20 * {T[1][M + 12]:.4f} + 0.80 * {T[1][M - 3]:.4f} = {p_beta_initial:.4f}")
    else:
        print("Option Beta: Not affordable.")
    
    print("\n--- Final Result ---")
    optimal_choice = "Alpha" if p_alpha_initial > p_beta_initial else "Beta"
    if p_alpha_initial == p_beta_initial and p_alpha_initial == 0.0:
        optimal_choice = "Neither (probability is 0)"
    elif p_alpha_initial == p_beta_initial:
        optimal_choice = "Alpha or Beta"

    print(f"Optimal first move is {optimal_choice}.")
    print(f"The maximum probability of success is: {final_probability}")

# --- Example Execution ---
# You can change these values to test different scenarios.
initial_investment_M = 25
number_of_trades_n = 3

solve_trading_problem(initial_investment_M, number_of_trades_n)
import sys

# It is recommended to increase Python's recursion limit for deep DP states, though our iterative approach avoids this.
# sys.setrecursionlimit(2000)

def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    by executing exactly n trades, using dynamic programming.
    """
    if M <= 0 or n <= 0:
        print("Initial investment and number of trades must be positive.")
        return

    print(f"Initial Investment (M): £{M}")
    print(f"Number of Trades (n): {n}")
    print(f"Target Investment: £{2 * M}\n")

    # The maximum possible money we can have is M + 12*n (n successful Beta trades).
    max_money = M + 12 * n

    # DP table: T[l][m] is the max probability of success with l trades remaining and m money.
    # l ranges from 0 to n.
    # m ranges from 0 to max_money.
    T = [[0.0 for _ in range(max_money + 1)] for _ in range(n + 1)]

    # Base case: l = 0 trades remaining.
    # Success is 1.0 if we have exactly 2*M, otherwise 0.
    if 2 * M <= max_money:
        T[0][2 * M] = 1.0

    # Fill the DP table using backward induction
    # l: number of trades remaining
    for l in range(1, n + 1):
        # m: current amount of money
        for m in range(max_money + 1):
            prob_alpha = 0.0
            # Strategy Alpha is possible if m >= 1
            if m >= 1:
                # The money after the trade will be m-1+2=m+1 or m-1+0=m-1
                prob_alpha = 0.6 * T[l - 1][m + 1] + 0.4 * T[l - 1][m - 1]

            prob_beta = 0.0
            # Strategy Beta is possible if m >= 3
            if m >= 3:
                # Check for boundary conditions, although with our max_money range it should be safe.
                m_beta_success = m + 12
                m_beta_failure = m - 3
                
                # Retrieve probabilities from the previous state (l-1)
                p_success = T[l-1][m_beta_success] if m_beta_success <= max_money else 0.0
                p_failure = T[l-1][m_beta_failure] if m_beta_failure >= 0 else 0.0

                prob_beta = 0.2 * p_success + 0.8 * p_failure

            T[l][m] = max(prob_alpha, prob_beta)

    final_probability = T[n][M]
    
    print("--- Optimal Strategy Analysis ---")
    if final_probability == 0.0:
        print("It's impossible to double the investment under the given constraints.")
        return

    print(f"The maximum probability of reaching £{2*M} is: {final_probability:.6f}\n")
    print("Analysis for the first trade:")
    
    # Show the calculation for the final state T[n][M]
    # Re-calculate the first step to show the breakdown
    
    # Strategy Alpha from initial state (M, n)
    prob_alpha_final = 0.0
    prob_alpha_succ_term, prob_alpha_fail_term = 0, 0
    if M >= 1:
        prob_alpha_succ_term = T[n - 1][M + 1]
        prob_alpha_fail_term = T[n - 1][M - 1]
        prob_alpha_final = 0.6 * prob_alpha_succ_term + 0.4 * prob_alpha_fail_term
        print(f"If we choose Strategy Alpha (cost £1):")
        print(f"  P(success) = 0.60 * P(at £{M+1}, {n-1} trades left) + 0.40 * P(at £{M-1}, {n-1} trades left)")
        print(f"  P(success) = 0.60 * {prob_alpha_succ_term:.4f} + 0.40 * {prob_alpha_fail_term:.4f} = {prob_alpha_final:.6f}")
    else:
        print("Strategy Alpha is not affordable.")

    # Strategy Beta from initial state (M, n)
    prob_beta_final = 0.0
    prob_beta_succ_term, prob_beta_fail_term = 0, 0
    if M >= 3:
        prob_beta_succ_term = T[n - 1][M + 12]
        prob_beta_fail_term = T[n - 1][M - 3]
        prob_beta_final = 0.2 * prob_beta_succ_term + 0.8 * prob_beta_fail_term
        print(f"If we choose Strategy Beta (cost £3):")
        print(f"  P(success) = 0.20 * P(at £{M+12}, {n-1} trades left) + 0.80 * P(at £{M-3}, {n-1} trades left)")
        print(f"  P(success) = 0.20 * {prob_beta_succ_term:.4f} + 0.80 * {prob_beta_fail_term:.4f} = {prob_beta_final:.6f}")
    else:
        print("Strategy Beta is not affordable.")
    
    print("\n--- Conclusion ---")
    if prob_alpha_final > prob_beta_final:
        print("The optimal first move is Strategy Alpha.")
    elif prob_beta_final > prob_alpha_final:
        print("The optimal first move is Strategy Beta.")
    else:
        # This handles cases where one is not affordable or probabilities are equal
        if prob_alpha_final == 0 and prob_beta_final == 0:
            print("No profitable path exists from the start.")
        else:
            print("Both strategies yield the same probability. Either can be chosen.")


if __name__ == '__main__':
    # --- Parameters for the simulation ---
    # You can change these values to test different scenarios.
    INITIAL_MONEY_M = 25
    NUMBER_OF_TRADES_n = 5
    
    solve_trading_problem(INITIAL_MONEY_M, NUMBER_OF_TRADES_n)

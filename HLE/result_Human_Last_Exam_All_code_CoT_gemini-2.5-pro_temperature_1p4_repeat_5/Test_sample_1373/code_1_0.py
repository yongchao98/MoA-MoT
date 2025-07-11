import sys

def solve_trading_problem():
    """
    This function designs and solves the trading problem using dynamic programming.
    """
    # Problem parameters (example values)
    M = 11
    n = 2

    print(f"Calculating the optimal probability for an initial investment M = {M} and n = {n} trades.")
    print(f"The target is to reach exactly 2*M = {2*M} GBP.\n")

    # Step 2: DP Table Setup
    max_capital = M + 12 * n
    target_capital = 2 * M

    # T[l][m]: probability of success from state (l trades done, m capital)
    T = [[0.0 for _ in range(max_capital + 1)] for _ in range(n + 1)]
    # strategy[l][m]: stores the best strategy choice ('A' for Alpha, 'B' for Beta)
    strategy = [['' for _ in range(max_capital + 1)] for _ in range(n + 1)]

    # Step 3: Base Case
    if target_capital <= max_capital:
        T[n][target_capital] = 1.0
    else:
        # Target is unreachable, probability is 0.
        print(f"The target capital {target_capital} is impossible to reach from M={M} in n={n} trades.")
        print("Maximum possible capital is:", max_capital)
        print("The probability of success is 0.")
        return

    # Step 4: Recurrence Relation (Iterate backwards)
    for l in range(n - 1, -1, -1):
        for m in range(max_capital + 1):
            prob_alpha = -1.0 # Use -1 to indicate not possible/calculated
            if m >= 1: # Check if we can afford Alpha
                # Capital after successful Alpha: m - 1 (fee) + 2 (return) = m + 1
                # Capital after failed Alpha:    m - 1 (fee) + 0 (return) = m - 1
                p_success = T[l+1][m+1] if m + 1 <= max_capital else 0.0
                p_fail = T[l+1][m-1] if m - 1 >= 0 else 0.0
                prob_alpha = 0.60 * p_success + 0.40 * p_fail

            prob_beta = -1.0
            if m >= 3: # Check if we can afford Beta
                # Capital after successful Beta: m - 3 (fee) + 15 (return) = m + 12
                # Capital after failed Beta:    m - 3 (fee) + 0 (return) = m - 3
                p_success = T[l+1][m+12] if m + 12 <= max_capital else 0.0
                p_fail = T[l+1][m-3] if m - 3 >= 0 else 0.0
                prob_beta = 0.20 * p_success + 0.80 * p_fail

            # Choose the optimal strategy
            if prob_alpha >= prob_beta:
                T[l][m] = prob_alpha
                strategy[l][m] = 'A'
            elif prob_beta > prob_alpha:
                T[l][m] = prob_beta
                strategy[l][m] = 'B'
            else: # both are negative (not possible)
                T[l][m] = 0

    # Step 5: Final Answer
    final_prob = T[0][M]
    initial_strategy = strategy[0][M]

    print("--- Results ---")
    print(f"The maximum probability of success is: {final_prob:.4f}")

    # Outputting the numbers in the final equation as requested
    print("\nThis result is derived from the first optimal trade:")
    if final_prob > 0 and initial_strategy:
        if initial_strategy == 'A':
            p_success_val = T[1][M+1] if M + 1 <= max_capital else 0.0
            p_fail_val = T[1][M-1] if M - 1 >= 0 else 0.0
            print(f"Optimal first move is Strategy Alpha from capital M={M}.")
            print(f"Prob(Success) = 0.60 * Prob(Success | l=1, m={M+1}) + 0.40 * Prob(Success | l=1, m={M-1})")
            print(f"Prob(Success) = 0.60 * {p_success_val:.4f} + 0.40 * {p_fail_val:.4f}")
            print(f"Prob(Success) = {0.60 * p_success_val + 0.40 * p_fail_val:.4f}")
        elif initial_strategy == 'B':
            p_success_val = T[1][M+12] if M + 12 <= max_capital else 0.0
            p_fail_val = T[1][M-3] if M - 3 >= 0 else 0.0
            print(f"Optimal first move is Strategy Beta from capital M={M}.")
            print(f"Prob(Success) = 0.20 * Prob(Success | l=1, m={M+12}) + 0.80 * Prob(Success | l=1, m={M-3})")
            print(f"Prob(Success) = 0.20 * {p_success_val:.4f} + 0.80 * {p_fail_val:.4f}")
            print(f"Prob(Success) = {0.20 * p_success_val + 0.80 * p_fail_val:.4f}")
    else:
        print("No profitable strategy sequence exists from the start.")

# Execute the solution
solve_trading_problem()
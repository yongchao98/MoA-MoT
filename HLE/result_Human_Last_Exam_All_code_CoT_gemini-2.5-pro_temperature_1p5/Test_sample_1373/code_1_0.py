import numpy as np

def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    to 2M in exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The exact number of trades to be executed.

    Returns:
        None. Prints the result and the final calculation.
    """
    if M <= 0 or n < 0:
        print("Initial investment and number of trades must be positive.")
        return

    # Maximum possible money: starting M + n successful Beta trades (profit £12 each)
    max_money = M + 12 * n
    target_money = 2 * M

    # DP table: T[l][m] = max probability with l trades left and m money
    # Dimensions: (n+1) rows for trades, (max_money+1) columns for money
    T = np.zeros((n + 1, max_money + 1))

    # Base case: l = 0 trades remaining
    if target_money <= max_money:
        T[0][target_money] = 1.0

    # Fill the DP table bottom-up for l = 1 to n trades remaining
    for l in range(1, n + 1):
        for m in range(max_money + 1):
            # Option 1: Strategy Alpha (Cost £1)
            p_alpha = 0.0
            if m >= 1:
                # States after trade: m+1 (success) and m-1 (failure)
                prob_success_alpha = T[l - 1][m + 1] if m + 1 <= max_money else 0.0
                prob_failure_alpha = T[l - 1][m - 1] if m - 1 >= 0 else 0.0
                p_alpha = 0.6 * prob_success_alpha + 0.4 * prob_failure_alpha

            # Option 2: Strategy Beta (Cost £3)
            p_beta = 0.0
            if m >= 3:
                # States after trade: m+12 (success) and m-3 (failure)
                prob_success_beta = T[l - 1][m + 12] if m + 12 <= max_money else 0.0
                prob_failure_beta = T[l - 1][m - 3] if m - 3 >= 0 else 0.0
                p_beta = 0.2 * prob_success_beta + 0.8 * prob_failure_beta
            
            # Determine T[l][m] based on the optimal choice
            if m < 1:
                T[l][m] = 0  # Cannot afford any trade
            elif m < 3:
                T[l][m] = p_alpha  # Only Alpha is possible
            else:
                T[l][m] = max(p_alpha, p_beta) # Choose the best strategy

    final_probability = T[n][M]

    print(f"Initial Investment (M): £{M}")
    print(f"Number of Trades (n): {n}")
    print(f"Target Investment (2M): £{target_money}")
    print("-" * 30)
    print(f"Maximum probability of success: {final_probability:.6f}\n")

    # Output the final calculation step as requested
    print(f"--- Breakdown of Final Calculation for T[{n}][{M}] ---")
    
    # Recalculate P_alpha and P_beta for the final state to print the equation
    p_alpha_final = 0.0
    val_alpha_succ = 0.0
    val_alpha_fail = 0.0
    if M >= 1:
        val_alpha_succ = T[n - 1][M + 1] if M + 1 <= max_money else 0.0
        val_alpha_fail = T[n - 1][M - 1] if M - 1 >= 0 else 0.0
        p_alpha_final = 0.6 * val_alpha_succ + 0.4 * val_alpha_fail
        print("Strategy Alpha:")
        print(f"  P = 0.60 * T[{n-1}][{M+1}] + 0.40 * T[{n-1}][{M-1}]")
        print(f"  P = 0.60 * {val_alpha_succ:.4f} + 0.40 * {val_alpha_fail:.4f} = {p_alpha_final:.6f}")
    else:
        print("Strategy Alpha: Not possible (not enough funds).")

    p_beta_final = 0.0
    val_beta_succ = 0.0
    val_beta_fail = 0.0
    if M >= 3:
        val_beta_succ = T[n - 1][M + 12] if M + 12 <= max_money else 0.0
        val_beta_fail = T[n - 1][M - 3] if M - 3 >= 0 else 0.0
        p_beta_final = 0.2 * val_beta_succ + 0.8 * val_beta_fail
        print("Strategy Beta:")
        print(f"  P = 0.20 * T[{n-1}][{M+12}] + 0.80 * T[{n-1}][{M-3}]")
        print(f"  P = 0.20 * {val_beta_succ:.4f} + 0.80 * {val_beta_fail:.4f} = {p_beta_final:.6f}")
    else:
        print("Strategy Beta: Not possible (not enough funds).")

    if final_probability > 0:
      if p_alpha_final > p_beta_final:
          print("\nOptimal choice for the first trade is Alpha.")
      elif p_beta_final > p_alpha_final:
          print("\nOptimal choice for the first trade is Beta.")
      else:
          print("\nBoth strategies yield the same probability.")


if __name__ == '__main__':
    # Example values for M and n
    initial_investment_M = 25
    num_trades_n = 10
    solve_trading_probability(initial_investment_M, num_trades_n)
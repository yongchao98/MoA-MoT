import collections

def solve_trading_probability(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    after exactly n trades using dynamic programming.

    Args:
        M (int): The initial investment in GBP.
        n (int): The number of trades to be executed.
    """
    if M <= 0 or n <= 0:
        print("Initial investment and number of trades must be positive.")
        return

    target_money = 2 * M

    # dp[i] is a dictionary mapping money -> max_probability for a state with i trades remaining.
    # We use defaultdict(float) to handle non-existent keys, which have 0 probability.
    dp = [collections.defaultdict(float) for _ in range(n + 1)]

    # Base Case: i = 0 trades remaining.
    # Probability of success is 1.0 only if we have exactly the target money.
    dp[0][target_money] = 1.0

    # Iterate backwards from i = 1 to n trades remaining.
    for i in range(1, n + 1):
        # We only need to calculate probabilities for money amounts 'm'
        # from which it's possible to reach the money amounts 'm_prev' that had
        # a non-zero success probability in the previous step (i-1).
        money_to_check = set()
        for m_prev in dp[i - 1]:
            # Potential previous states 'm' for Strategy Alpha
            money_to_check.add(m_prev - 1)  # successful alpha came from m_prev - 1
            money_to_check.add(m_prev + 1)  # failed alpha came from m_prev + 1
            # Potential previous states 'm' for Strategy Beta
            money_to_check.add(m_prev - 12) # successful beta came from m_prev - 12
            money_to_check.add(m_prev + 3)  # failed beta came from m_prev + 3

        for m in money_to_check:
            # We cannot have negative money or pay fees we can't afford.
            if m < 0:
                continue

            # --- Calculate probability if we choose Strategy Alpha ---
            prob_alpha = 0.0
            if m >= 1:
                prob_success_from_win = dp[i - 1][m + 1]
                prob_success_from_loss = dp[i - 1][m - 1]
                prob_alpha = (0.6 * prob_success_from_win) + (0.4 * prob_success_from_loss)

            # --- Calculate probability if we choose Strategy Beta ---
            prob_beta = 0.0
            if m >= 3:
                prob_success_from_win = dp[i - 1][m + 12]
                prob_success_from_loss = dp[i - 1][m - 3]
                prob_beta = (0.2 * prob_success_from_win) + (0.8 * prob_success_from_loss)

            # The optimal strategy for state (i, m) maximizes the probability.
            dp[i][m] = max(prob_alpha, prob_beta)

    # The final answer is the probability starting with M pounds and n trades.
    final_probability = dp[n][M]

    # --- Output the results and the breakdown for the first decision ---
    print(f"Initial Investment M: {M} GBP")
    print(f"Number of Trades n: {n}")
    print(f"Target Investment: {target_money} GBP")
    print("-" * 40)
    print(f"The maximum probability of reaching exactly {target_money} GBP is: {final_probability:.7f}")
    print("-" * 40)
    
    # To fulfill the prompt's request, let's show the numbers in the final calculation for the first trade
    print("Calculation for the Optimal First Trade:")
    print(f"State: {n} trades remaining, {M} GBP available.")
    
    prob_alpha_final = 0.0
    if M >= 1:
        prob_s_a_w = dp[n-1][M+1]
        prob_s_a_l = dp[n-1][M-1]
        prob_alpha_final = (0.6 * prob_s_a_w) + (0.4 * prob_s_a_l)
        print(f"\nStrategy Alpha (Cost: 1 GBP):")
        print(f"  P(success) = 0.60 * P(success | {n-1} trades, {M+1} GBP) + 0.40 * P(success | {n-1} trades, {M-1} GBP)")
        print(f"  P(success) = 0.60 * {prob_s_a_w:.5f} + 0.40 * {prob_s_a_l:.5f}")
        print(f"  Resulting Probability: {prob_alpha_final:.7f}")
    else:
        print("\nStrategy Alpha: Not affordable.")

    prob_beta_final = 0.0
    if M >= 3:
        prob_s_b_w = dp[n-1][M+12]
        prob_s_b_l = dp[n-1][M-3]
        prob_beta_final = (0.2 * prob_s_b_w) + (0.8 * prob_s_b_l)
        print(f"\nStrategy Beta (Cost: 3 GBP):")
        print(f"  P(success) = 0.20 * P(success | {n-1} trades, {M+12} GBP) + 0.80 * P(success | {n-1} trades, {M-3} GBP)")
        print(f"  P(success) = 0.20 * {prob_s_b_w:.5f} + 0.80 * {prob_s_b_l:.5f}")
        print(f"  Resulting Probability: {prob_beta_final:.7f}")
    else:
        print("\nStrategy Beta: Not affordable.")

    print("-" * 40)
    if final_probability > 0:
      if prob_alpha_final > prob_beta_final:
          print("Optimal first move is Strategy Alpha.")
      elif prob_beta_final > prob_alpha_final:
          print("Optimal first move is Strategy Beta.")
      else:
          print("Optimal first move can be either Alpha or Beta (they yield the same probability).")
    else:
        print("There is no strategy that results in a non-zero probability of success.")


# Example Usage:
# You can change these values to test different scenarios.
initial_investment_M = 10
number_of_trades_n = 5
solve_trading_probability(initial_investment_M, number_of_trades_n)
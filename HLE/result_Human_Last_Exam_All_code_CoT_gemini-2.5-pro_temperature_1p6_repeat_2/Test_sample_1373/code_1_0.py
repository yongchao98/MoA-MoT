import sys

# It's possible to hit the recursion limit for large n if using a recursive solution.
# For the iterative approach here, this is not strictly necessary but good practice.
sys.setrecursionlimit(2000)

def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    to 2M in exactly n trades using dynamic programming.
    """
    target_money = 2 * M

    # dp[k][m] will store a tuple: (max_probability, optimal_choice)
    # k: trades remaining
    # m: money
    # Using dictionaries for sparse storage of money states.
    dp = {}

    # Base Case: k = 0 trades remaining
    dp[0] = {target_money: (1.0, 'Success')}

    # Iterate from k = 1 to n trades remaining
    for k in range(1, n + 1):
        dp[k] = {}
        # We only need to compute probabilities for money amounts from which
        # it's possible to reach the states with non-zero probability in the previous step (k-1).
        money_to_compute = set()
        if k - 1 in dp:
            for m_prev in dp[k - 1]:
                # Inverse transitions from m_prev to find potential m values for this step
                # Alpha: m -> m_prev => m+1=m_prev or m-1=m_prev
                money_to_compute.add(m_prev - 1)
                money_to_compute.add(m_prev + 1)
                # Beta: m -> m_prev => m+12=m_prev or m-3=m_prev
                money_to_compute.add(m_prev - 12)
                money_to_compute.add(m_prev + 3)

        for m in sorted(list(money_to_compute)):
            # We cannot have negative money
            if m < 0:
                continue

            # --- Calculate probability from Strategy Alpha ---
            prob_alpha = -1.0 # Use -1 to indicate not a valid option
            if m >= 1:
                # Get the probability of success from the resulting states in the next step (k-1)
                prob_win_alpha = dp.get(k - 1, {}).get(m + 1, (0.0, None))[0]
                prob_lose_alpha = dp.get(k - 1, {}).get(m - 1, (0.0, None))[0]
                prob_alpha = 0.6 * prob_win_alpha + 0.4 * prob_lose_alpha

            # --- Calculate probability from Strategy Beta ---
            prob_beta = -1.0 # Use -1 to indicate not a valid option
            if m >= 3:
                # Get the probability of success from the resulting states in the next step (k-1)
                prob_win_beta = dp.get(k - 1, {}).get(m + 12, (0.0, None))[0]
                prob_lose_beta = dp.get(k - 1, {}).get(m - 3, (0.0, None))[0]
                prob_beta = 0.2 * prob_win_beta + 0.8 * prob_lose_beta

            # --- Determine the optimal strategy and store the result ---
            if prob_alpha > prob_beta:
                dp[k][m] = (prob_alpha, 'Alpha')
            elif prob_beta >= 0 and prob_beta >= prob_alpha:
                dp[k][m] = (prob_beta, 'Beta')
            # If neither strategy is chosen or possible, it won't be added to dp[k]

    # --- Final Result and Equation ---
    final_prob, final_choice = dp.get(n, {}).get(M, (0.0, None))

    print(f"Initial investment M = £{M}, Trades n = {n}, Target = £{2*M}\n")
    print(f"Let T(m, k) be the optimal probability of success with £m and k trades remaining.")
    print(f"\nThe maximum probability of success is: T({M}, {n}) = {final_prob:.7f}\n")

    if final_prob > 0:
        print("This result is derived from the following top-level calculation:")
        if final_choice == 'Alpha':
            m_win, m_lose = M + 1, M - 1
            p_win = dp.get(n - 1, {}).get(m_win, (0.0, None))[0]
            p_lose = dp.get(n - 1, {}).get(m_lose, (0.0, None))[0]
            
            print(f"The optimal choice for the first trade is Strategy Alpha.")
            print(f"T({M}, {n}) = 0.6 * T({m_win}, {n-1}) + 0.4 * T({m_lose}, {n-1})")
            print(f"T({M}, {n}) = 0.6 * {p_win:.7f} + 0.4 * {p_lose:.7f}")
            print(f"T({M}, {n}) = {0.6 * p_win:.7f} + {0.4 * p_lose:.7f} = {final_prob:.7f}")
        
        elif final_choice == 'Beta':
            m_win, m_lose = M + 12, M - 3
            p_win = dp.get(n - 1, {}).get(m_win, (0.0, None))[0]
            p_lose = dp.get(n - 1, {}).get(m_lose, (0.0, None))[0]
            
            print(f"The optimal choice for the first trade is Strategy Beta.")
            print(f"T({M}, {n}) = 0.2 * T({m_win}, {n-1}) + 0.8 * T({m_lose}, {n-1})")
            print(f"T({M}, {n}) = 0.2 * {p_win:.7f} + 0.8 * {p_lose:.7f}")
            print(f"T({M}, {n}) = {0.2 * p_win:.7f} + {0.8 * p_lose:.7f} = {final_prob:.7f}")
    else:
        print("It is impossible to guarantee reaching the target under any strategy.")


if __name__ == '__main__':
    # --- Example Usage ---
    # You can change these values to test different scenarios
    initial_investment_M = 25
    number_of_trades_n = 5
    
    solve_trading_problem(initial_investment_M, number_of_trades_n)

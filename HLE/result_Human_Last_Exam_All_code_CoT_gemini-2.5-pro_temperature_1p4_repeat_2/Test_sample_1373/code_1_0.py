import sys

# It's good practice to increase recursion limit for deep DP, although not strictly necessary for this iterative solution.
sys.setrecursionlimit(2000)

def solve_trading_problem(M, n):
    """
    Calculates the maximum probability of doubling an initial investment M
    in exactly n trades using dynamic programming.
    """
    if M <= 0 or n <= 0:
        print("Initial investment and number of trades must be positive.")
        return

    TARGET = 2 * M
    # The maximum possible amount of money we can have.
    # Start with M, max gain per trade is +12 from Beta.
    MAX_MONEY = M + 12 * n

    # dp[i][j]: max probability of success with i trades left and j pounds.
    # Initialize with 0.0
    dp = [[0.0 for _ in range(MAX_MONEY + 1)] for _ in range(n + 1)]
    
    # choices[i][j]: optimal choice at state (i, j)
    choices = [['none' for _ in range(MAX_MONEY + 1)] for _ in range(n + 1)]

    # Base case: 0 trades left. Success only if money is exactly TARGET.
    if TARGET <= MAX_MONEY:
        dp[0][TARGET] = 1.0

    # Fill DP table bottom-up (or trades-left up)
    for i in range(1, n + 1):  # i = number of trades left
        for j in range(MAX_MONEY + 1):  # j = current money
            
            p_alpha = 0.0
            p_beta = 0.0

            # --- Calculate probability from choosing Strategy Alpha ---
            # Requires £1 fee
            if j >= 1:
                # 60% chance of success (net +£1)
                prob_alpha_s = 0.0
                if j + 1 <= MAX_MONEY:
                    prob_alpha_s = dp[i-1][j+1]
                
                # 40% chance of failure (net -£1)
                prob_alpha_f = 0.0
                if j - 1 >= 0:
                    prob_alpha_f = dp[i-1][j-1]

                p_alpha = 0.60 * prob_alpha_s + 0.40 * prob_alpha_f
            
            # --- Calculate probability from choosing Strategy Beta ---
            # Requires £3 fee
            if j >= 3:
                # 20% chance of success (net +£12)
                prob_beta_s = 0.0
                if j + 12 <= MAX_MONEY:
                    prob_beta_s = dp[i-1][j+12]

                # 80% chance of failure (net -£3)
                prob_beta_f = 0.0
                if j - 3 >= 0:
                    prob_beta_f = dp[i-1][j-3]

                p_beta = 0.20 * prob_beta_s + 0.80 * prob_beta_f

            # --- Decide optimal strategy ---
            if p_alpha > p_beta:
                dp[i][j] = p_alpha
                choices[i][j] = 'Alpha'
            else:
                # If p_beta is greater or they are equal, we prefer Beta
                # (or it doesn't matter if both are 0).
                dp[i][j] = p_beta
                if p_beta > 0:
                    choices[i][j] = 'Beta'

    # --- Final result and explanation ---
    final_prob = dp[n][M]
    initial_choice = choices[n][M]

    print(f"Initial State: Investment M = £{M}, Trades n = {n}")
    print(f"Target Amount: £{TARGET}")
    print("-" * 30)

    if final_prob == 0.0:
        print("It's impossible to guarantee reaching the target amount under any strategy.")
    else:
        print(f"The optimal strategy for the first trade is: Strategy {initial_choice}")
        print("The calculation for the initial state is:")
        
        if initial_choice == 'Alpha':
            # Note: fee=1, return=2 -> net=+1
            # fee=1, return=0 -> net=-1
            val1 = dp[n-1][M+1]
            val2 = dp[n-1][M-1]
            print(f"P(success) = 0.60 * P(success | {n-1} trades left, £{M+1}) + 0.40 * P(success | {n-1} trades left, £{M-1})")
            print(f"P(success) = 0.60 * {val1:.4f} + 0.40 * {val2:.4f}")
        elif initial_choice == 'Beta':
            # Note: fee=3, return=15 profit -> net=+12
            # fee=3, return=0 profit -> net=-3
            val1 = dp[n-1][M+12]
            val2 = dp[n-1][M-3]
            print(f"P(success) = 0.20 * P(success | {n-1} trades left, £{M+12}) + 0.80 * P(success | {n-1} trades left, £{M-3})")
            print(f"P(success) = 0.20 * {val1:.4f} + 0.80 * {val2:.4f}")

    print("-" * 30)
    print(f"Maximum probability of reaching exactly £{TARGET}: {final_prob:.6f}")


# Example Usage: Let's use an example where a non-zero probability exists.
# A net profit of M=14 can be achieved in n=3 trades (e.g., Beta success, Alpha success, Alpha success).
# Profit = +12 + 1 + 1 = 14.
INITIAL_MONEY = 14
NUM_TRADES = 3
solve_trading_problem(INITIAL_MONEY, NUM_TRADES)

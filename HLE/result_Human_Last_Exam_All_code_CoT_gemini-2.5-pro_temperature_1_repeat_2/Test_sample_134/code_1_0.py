import numpy as np

def solve_coin_game():
    """
    Solves the coin line game for a given arrangement using dynamic programming.
    """
    # Define the number of coins
    num_1_euro = 136
    num_2_euro = 87
    total_coins = num_1_euro + num_2_euro

    # Create a "worst-case" arrangement for Player 1.
    # We place the 2-euro coins in the middle. To give P2 the advantage on this
    # subgame, the number of 1-euro coins picked before reaching the 2s should be odd.
    # Let's put 69 1-euro coins on the left, and 136-69=67 on the right.
    # 69 coins will be picked (35 by P1, 34 by P2). Turn 70 is P2's turn to start on the 2s.
    coins = [1] * 69 + [2] * num_2_euro + [1] * 67
    
    n = len(coins)
    
    # dp[i][j] stores the maximum value the current player can get from the sub-array coins[i:j+1]
    dp = np.zeros((n, n))
    
    # Pre-calculate prefix sums to get S(i,j) in O(1)
    prefix_sum = [0] * (n + 1)
    for i in range(n):
        prefix_sum[i+1] = prefix_sum[i] + coins[i]

    def get_sum(i, j):
        if i > j:
            return 0
        return prefix_sum[j+1] - prefix_sum[i]

    # Fill the DP table
    for length in range(1, n + 1):
        for i in range(n - length + 1):
            j = i + length - 1
            if i == j:
                dp[i][j] = coins[i]
                continue
            
            # If player takes coins[i]:
            # The other player will play on coins[i+1...j] and get dp[i+1][j].
            # The current player's total score will be coins[i] + (sum of remaining - other player's score)
            val1 = coins[i] + get_sum(i + 1, j) - dp[i + 1][j]
            
            # If player takes coins[j]:
            val2 = coins[j] + get_sum(i, j - 1) - dp[i][j - 1]
            
            dp[i][j] = max(val1, val2)

    p1_score = dp[0][n-1]
    total_value = get_sum(0, n - 1)
    p2_score = total_value - p1_score

    print(f"Game setup: {num_1_euro} 1-euro coins and {num_2_euro} 2-euro coins.")
    print(f"Total value on the line: {total_value} euros.")
    print("Analyzing a specific arrangement designed to be challenging for Player 1...")
    print(f"Player 1's optimal score: {int(p1_score)}")
    print(f"Player 2's optimal score: {int(p2_score)}")

    if p1_score > p2_score:
        print("Result: Player 1 wins.")
    elif p2_score > p1_score:
        print("Result: Player 2 wins.")
    else:
        print("Result: It's a draw.")

solve_coin_game()
<<<A>>>
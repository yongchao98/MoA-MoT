import random
import sys

# Increase recursion limit for the solver
sys.setrecursionlimit(300)

def solve_game():
    """
    Analyzes the coin game, demonstrating Player 1's advantage.
    """
    # 1. Game Setup
    num_1_euro = 136
    num_2_euro = 87
    coins = [1] * num_1_euro + [2] * num_2_euro
    random.shuffle(coins)
    
    n = len(coins)
    total_value = sum(coins)
    
    print(f"A random line of {n} coins has been arranged.")
    print(f"Total value of all coins is: {total_value} euros.")
    print("-" * 30)

    # Memoization cache for the solver
    memo = {}
    # Pre-compute prefix sums for efficient calculation of sub-array sums
    prefix_sum = [0] * (n + 1)
    for i in range(n):
        prefix_sum[i+1] = prefix_sum[i] + coins[i]

    def get_sum(i, j):
        """Calculates sum of coins from index i to j using prefix sums."""
        if i > j:
            return 0
        return prefix_sum[j+1] - prefix_sum[i]

    def find_max_score(i, j):
        """
        Recursively finds the maximum score the current player can get
        from the subarray of coins from index i to j.
        Uses memoization to avoid re-computing results.
        """
        if i > j:
            return 0
        if i == j:
            return coins[i]
        
        # Check cache
        if (i, j) in memo:
            return memo[(i, j)]

        # If the current player takes the coin at index i:
        # The other player will play on sub-array (i+1, j) and get find_max_score(i+1, j).
        # The current player's total score will be the coin they took plus what's left over.
        score_if_take_left = coins[i] + get_sum(i+1, j) - find_max_score(i+1, j)

        # If the current player takes the coin at index j:
        score_if_take_right = coins[j] + get_sum(i, j-1) - find_max_score(i, j-1)
        
        # The player chooses the move that maximizes their score.
        result = max(score_if_take_left, score_if_take_right)
        memo[(i, j)] = result
        return result

    # 3. Player 1's Analysis
    print("Player 1 analyzes the two possible first moves:")

    # --- Option A: Player 1 takes the leftmost coin ---
    p1_takes_coin_A = coins[0]
    # P2 plays on the rest of the line and gets the optimal score for that sub-game.
    p2_score_A = find_max_score(1, n - 1)
    # P1 gets what's left from the sub-game.
    p1_score_A = p1_takes_coin_A + (get_sum(1, n - 1) - p2_score_A)
    
    print(f"\nOption A: P1 takes the leftmost coin ({p1_takes_coin_A} euro).")
    print(f"  - P2 will then play optimally and get {p2_score_A} euros.")
    print(f"  - P1's final score would be: {p1_takes_coin_A} + ({get_sum(1, n - 1)} - {p2_score_A}) = {p1_score_A} euros.")

    # --- Option B: Player 1 takes the rightmost coin ---
    p1_takes_coin_B = coins[n-1]
    # P2 plays on the rest of the line and gets the optimal score for that sub-game.
    p2_score_B = find_max_score(0, n - 2)
    # P1 gets what's left from the sub-game.
    p1_score_B = p1_takes_coin_B + (get_sum(0, n - 2) - p2_score_B)

    print(f"\nOption B: P1 takes the rightmost coin ({p1_takes_coin_B} euro).")
    print(f"  - P2 will then play optimally and get {p2_score_B} euros.")
    print(f"  - P1's final score would be: {p1_takes_coin_B} + ({get_sum(0, n - 2)} - {p2_score_B}) = {p1_score_B} euros.")
    
    # 4. Player 1's Decision and Conclusion
    print("-" * 30)
    print("Decision:")

    if p1_score_A >= p1_score_B:
        print(f"Player 1 chooses Option A, guaranteeing a score of {p1_score_A}.")
        print(f"Player 2's score will be {p2_score_A}.")
        if p1_score_A > p2_score_A:
            print("Player 1 wins.")
        elif p2_score_A > p1_score_A:
            print("Player 2 wins.")
        else:
            print("The game is a tie.")
    else:
        print(f"Player 1 chooses Option B, guaranteeing a score of {p1_score_B}.")
        print(f"Player 2's score will be {p2_score_B}.")
        if p1_score_B > p2_score_B:
            print("Player 1 wins.")
        elif p2_score_B > p1_score_B:
            print("Player 2 wins.")
        else:
            print("The game is a tie.")

    print("\nBecause Player 1 can always perform this analysis and choose the better outcome for themselves, they have a clear strategic advantage.")
    print("Therefore, you should prefer to be the 1st player.")

if __name__ == '__main__':
    solve_game()
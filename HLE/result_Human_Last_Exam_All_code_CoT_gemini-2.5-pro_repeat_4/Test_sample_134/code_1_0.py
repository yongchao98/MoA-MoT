def solve_game_choice():
    """
    Analyzes the coin game to determine if it's better to be the 1st or 2nd player.
    """
    # Initial number of coins
    one_euro_coins = 136
    two_euro_coins = 87

    # --- Step 1: Calculate total coins ---
    total_coins = one_euro_coins + two_euro_coins
    print(f"Total number of coins: {one_euro_coins} + {two_euro_coins} = {total_coins}")

    # --- Step 2: Calculate coins per player ---
    # Since the total number of coins is odd, Player 1 gets one more turn.
    player_1_coins = (total_coins + 1) // 2
    player_2_coins = total_coins - player_1_coins
    print(f"Player 1 gets to pick {player_1_coins} coins.")
    print(f"Player 2 gets to pick {player_2_coins} coins.")

    # --- Step 3: Calculate total value ---
    total_value = (one_euro_coins * 1) + (two_euro_coins * 2)
    print(f"The total value on the table is: ({one_euro_coins} * 1) + ({two_euro_coins} * 2) = {total_value} euros.")
    
    # --- Step 4: Reasoning ---
    print("\nReasoning:")
    print("The total number of coins is odd (223). This means Player 1 will always pick one more coin than Player 2 (112 vs 111).")
    print("While a specific arrangement of coins might allow a skilled Player 2 to win by picking more valuable coins, Player 1 has the fundamental advantage of an extra turn.")
    print("Player 1 also makes the first move, allowing them to react to the initial layout and choose the move that leads to the best possible outcome for them, considering Player 2 will also play optimally.")
    print("Averaged over all possible random arrangements, the structural advantage lies with the first player due to the extra turn. Therefore, it is preferable to be the 1st player.")

solve_game_choice()
def solve_coin_game():
    """
    Analyzes the coin game to determine which player has the advantage.
    """
    # Define the number of each type of coin
    n_1_euro = 136
    n_2_euro = 87

    # Calculate the total number of coins
    total_coins = n_1_euro + n_2_euro

    # Since the total is odd, Player 1 gets one more coin than Player 2.
    player_1_coins = (total_coins + 1) // 2
    player_2_coins = total_coins // 2

    # Calculate the total value of all coins
    total_value = n_1_euro * 1 + n_2_euro * 2
    
    # The winning score is any amount greater than half the total value
    winning_threshold = total_value / 2

    print("--- Game Analysis ---")
    print(f"Total coins: {total_coins}")
    print(f"Player 1 gets to pick {player_1_coins} coins.")
    print(f"Player 2 gets to pick {player_2_coins} coins.")
    print(f"The total value of all coins is {n_1_euro} * 1 + {n_2_euro} * 2 = {total_value} euros.")
    print(f"To win, a player must collect a value greater than {winning_threshold} euros.")
    
    print("\n--- Conclusion ---")
    print("The key to this problem is that the total number of coins is odd (223).")
    print("This means Player 1, who goes first, will also take the last coin.")
    print(f"As a result, Player 1 collects {player_1_coins} coins while Player 2 collects only {player_2_coins}.")
    print("Because the choice of role is made before the coins are arranged randomly, we must base our decision on the game's structure, not a specific layout.")
    print("The structural advantage lies with Player 1, who gets an extra coin. This provides a fundamental advantage that cannot be reliably overcome by Player 2's strategy across all possible random arrangements.")
    print("\nTherefore, you should prefer to be the 1st player.")

solve_coin_game()
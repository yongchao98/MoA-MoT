def solve_coin_game():
    """
    Analyzes the coin game to determine if it's better to be the 1st or 2nd player.
    """
    
    # Step 1: Define the game parameters and calculate totals.
    num_1_euro = 136
    num_2_euro = 87
    
    print("Step 1: Analyzing the coin setup.")
    # Calculate and print the total number of coins
    total_coins = num_1_euro + num_2_euro
    print(f"The total number of coins is {num_1_euro} + {num_2_euro} = {total_coins}.")

    # Step 2: Determine the number of turns for each player.
    print("\nStep 2: Analyzing the number of turns per player.")
    print(f"The total number of coins ({total_coins}) is an odd number.")
    print("This means the first and second players will not pick the same number of coins.")
    
    # Calculate and print the turns for each player
    player1_turns = (total_coins + 1) // 2
    player2_turns = total_coins // 2
    print(f"Player 1 picks ({total_coins} + 1) / 2 = {player1_turns} coins.")
    print(f"Player 2 picks {total_coins} / 2 = {player2_turns} coins.")
    print("Player 1 gets to pick one more coin than Player 2.")

    # Step 3: Explain the winning strategy based on the game's structure.
    print("\nStep 3: Explaining the winning strategy for Player 1.")
    print("Because the number of coins is odd, Player 1 has a guaranteed winning strategy.")
    print("This strategy involves mentally pairing the coins from the ends inward:")
    print(f"(Coin at position 1, Coin at position {total_coins})")
    print(f"(Coin at position 2, Coin at position {total_coins - 1})")
    print("...")
    
    # Calculate the number of pairs and the middle coin position
    num_pairs = player2_turns
    middle_coin_pos = player1_turns
    print(f"This creates {num_pairs} pairs, with one coin left unpaired in the exact middle (at position {middle_coin_pos}).")
    
    print("\nOn every turn, Player 1 can choose which end to pick from. They can play to ensure that for each of the {} pairs, they get the coin with a value greater than or equal to the coin their opponent gets from that pair.".format(num_pairs))
    print(f"Furthermore, the middle coin at position {middle_coin_pos} will be the very last coin remaining, and it will be Player 1's turn to pick it.")
    
    # Step 4: Conclusion.
    print("\nStep 4: Conclusion.")
    print("By following this strategy, Player 1's total value will always be greater than Player 2's.")
    print("Therefore, it is always better to be the 1st player.")

# Run the analysis
solve_coin_game()
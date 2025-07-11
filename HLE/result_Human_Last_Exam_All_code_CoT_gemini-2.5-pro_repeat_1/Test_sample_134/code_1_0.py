def solve_coin_game():
    """
    Analyzes the coin game to determine if it's better to be the 1st or 2nd player.
    """
    # Step 1: Define the initial setup
    one_euro_coins = 136
    two_euro_coins = 87
    total_coins = one_euro_coins + two_euro_coins

    # Since the total number of coins is odd, Player 1 gets one more turn than Player 2.
    p1_turns = (total_coins + 1) // 2
    p2_turns = (total_coins - 1) // 2

    print(f"There are {one_euro_coins} 1-euro coins and {two_euro_coins} 2-euro coins, for a total of {total_coins} coins.")
    print(f"Player 1 will pick a total of {p1_turns} coins.")
    print(f"Player 2 will pick a total of {p2_turns} coins.")
    print("-" * 40)

    # Step 2 & 3: Analyze the score and the winning condition
    print("Let's analyze the final scores, V1 for Player 1 and V2 for Player 2.")
    print("Let N1_2 be the number of 2-euro coins Player 1 collects.")
    print("Let N2_2 be the number of 2-euro coins Player 2 collects.")
    print("\nThe score for a player is the sum of their coins' values.")
    print("We can express a player's score based on how many 2-euro coins they get.")
    
    print(f"\nFor Player 1: V1 = (1 * number of 1-euro coins) + (2 * N1_2)")
    print(f"Since Player 1 gets {p1_turns} coins, their number of 1-euro coins is ({p1_turns} - N1_2).")
    print(f"So, V1 = 1 * ({p1_turns} - N1_2) + 2 * N1_2 = {p1_turns} + N1_2")

    print(f"\nFor Player 2: V2 = (1 * number of 1-euro coins) + (2 * N2_2)")
    print(f"Since Player 2 gets {p2_turns} coins, their number of 1-euro coins is ({p2_turns} - N2_2).")
    print(f"So, V2 = 1 * ({p2_turns} - N2_2) + 2 * N2_2 = {p2_turns} + N2_2")

    print("\nTo see who wins, let's find the difference in scores, V1 - V2:")
    score_diff_constant = p1_turns - p2_turns
    # This is the final equation the user wants to see, with each number included.
    print(f"V1 - V2 = ({p1_turns} + N1_2) - ({p2_turns} + N2_2)")
    print(f"V1 - V2 = ({p1_turns} - {p2_turns}) + (N1_2 - N2_2)")
    print(f"The final, simplified equation is: V1 - V2 = {score_diff_constant} + N1_2 - N2_2")
    print("-" * 40)

    # Step 4: Assess the advantages and conclude
    print("This final equation reveals the core of the game:")
    print(f"Player 1 wins if V1 > V2, which means {score_diff_constant} + N1_2 - N2_2 > 0.")
    print(f"This simplifies to N1_2 > N2_2 - {score_diff_constant}, which means N1_2 >= N2_2 (since they are integers).")
    
    print("\nIn other words:")
    print("  - Player 1 wins if they collect at least as many 2-euro coins as Player 2.")
    print("  - A draw occurs if Player 2 gets exactly one more 2-euro coin than Player 1.")
    print("  - Player 2 wins only if they get two or more 2-euro coins than Player 1.")

    print("\nConclusion:")
    print("The game is essentially a race to collect the most 2-euro coins.")
    print("Player 1 has two structural advantages, regardless of how the coins are arranged:")
    print("1. They get the first move, allowing them to control the game from the start.")
    print(f"2. They get more turns ({p1_turns} vs {p2_turns}), giving them one extra opportunity to pick a coin.")
    print("\nThese advantages put Player 1 in a better position to win the 'race' for the 2-euro coins.")
    print("Therefore, you should prefer to be the 1st player.")

solve_coin_game()
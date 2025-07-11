def solve_coin_dilemma():
    """
    Analyzes the coin game to determine if it's better to be the 1st or 2nd player.
    """
    
    # 1. Define the game setup
    one_euro_coins = 136
    two_euro_coins = 87
    
    total_coins = one_euro_coins + two_euro_coins
    total_value = (one_euro_coins * 1) + (two_euro_coins * 2)
    
    print("--- Game Analysis ---")
    print(f"Total coins: {one_euro_coins} (1-euro) + {two_euro_coins} (2-euro) = {total_coins} coins")
    print(f"Total value: ({one_euro_coins} * 1) + ({two_euro_coins} * 2) = {total_value} euros")
    print(f"To win, a player needs a score greater than {total_value / 2} euros.")
    print("-" * 23)

    # 2. Analyze the turns
    # Since the total number of coins (223) is odd, the first player gets more turns.
    p1_coins_count = (total_coins + 1) // 2
    p2_coins_count = total_coins - p1_coins_count
    
    print("--- Turn Distribution ---")
    print(f"Because there are an odd number of coins ({total_coins}), players take a different number of coins:")
    print(f"Player 1 gets {p1_coins_count} coins.")
    print(f"Player 2 gets {p2_coins_count} coins.")
    print("Player 1 has the advantage of taking one extra coin.")
    print("-" * 23)

    # 3. Naive Expectation Analysis
    # While players play optimally, a simple expected value calculation reveals the underlying advantage.
    # It assumes the value of each coin picked is the average value of all coins in the line.
    average_coin_value = total_value / total_coins
    
    p1_expected_score = p1_coins_count * average_coin_value
    p2_expected_score = p2_coins_count * average_coin_value
    
    print("--- Expected Score Analysis ---")
    print("Although players play optimally, we can estimate the advantage by looking at the expected score.")
    print(f"The average value of a random coin is {total_value} / {total_coins} = {average_coin_value:.4f} euros.")
    print("\nIf each picked coin had this average value, the scores would be:")
    print(f"Player 1's Expected Score = {p1_coins_count} coins * {average_coin_value:.4f} = {p1_expected_score:.4f} euros")
    print(f"Player 2's Expected Score = {p2_coins_count} coins * {average_coin_value:.4f} = {p2_expected_score:.4f} euros")
    
    # 4. Conclusion
    print("\n--- Conclusion ---")
    print("Player 2 has strategic options to counter Player 1. For certain specific arrangements, Player 2 can win.")
    print("However, the problem states the coins are arranged at random.")
    print("Across random arrangements, the fundamental advantage of getting one more coin makes Player 1's position more favorable.")
    print("The expected score for Player 1 is higher and crosses the winning threshold of 155.")
    print("\nTherefore, it is preferable to be the 1st player.")

solve_coin_dilemma()
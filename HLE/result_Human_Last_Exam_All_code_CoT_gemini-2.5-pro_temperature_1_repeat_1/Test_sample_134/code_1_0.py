def solve_coin_game():
    """
    This function explains the reasoning behind choosing to be the 1st or 2nd player
    and calculates the win condition.
    """
    
    # Game setup
    one_euro_coins = 136
    two_euro_coins = 87
    
    # Total coins and turns
    total_coins = one_euro_coins + two_euro_coins
    
    # Number of coins each player gets
    p1_coins = (total_coins + 1) // 2
    p2_coins = (total_coins - 1) // 2
    
    # Total value of all coins
    total_value = (one_euro_coins * 1) + (two_euro_coins * 2)
    
    # Let n2_1 be the number of 2-euro coins Player 1 collects.
    # Player 1's score is V1 = (p1_coins - n2_1) * 1 + n2_1 * 2 = p1_coins + n2_1
    # Player 2's score is V2 = total_value - V1 = total_value - (p1_coins + n2_1)
    
    # We can express the win condition without knowing the exact value of n2_1.
    # V1 > V2  => p1_coins + n2_1 > total_value - p1_coins - n2_1
    #         => 2 * p1_coins + 2 * n2_1 > total_value
    #         => 2 * (p1_coins + n2_1) > total_value
    #
    # Let's find the threshold for n2_1 that results in a draw (V1 = V2).
    # A draw means V1 = V2 = total_value / 2
    # p1_coins + n2_1_draw = total_value / 2
    # n2_1_draw = (total_value / 2) - p1_coins
    
    n2_1_draw = (total_value / 2.0) - p1_coins
    
    # Player 1 wins if they get more 2-euro coins than this threshold.
    # n2_1 > n2_1_draw
    
    print("Step 1: Analyze the number of coins per player.")
    print(f"There are a total of {total_coins} coins, which is an odd number.")
    print(f"Player 1 will pick {p1_coins} coins.")
    print(f"Player 2 will pick {p2_coins} coins.")
    print("Player 1 has a fundamental advantage by getting one more coin.\n")

    print("Step 2: Determine the win condition based on the value.")
    print(f"The total value of all coins is ({one_euro_coins} * 1) + ({two_euro_coins} * 2) = {total_value} euros.")
    print("Let n2_1 be the number of 2-euro coins Player 1 collects.")
    print(f"Player 1's score is V1 = {p1_coins} + n2_1.")
    print(f"Player 2's score is V2 = {total_value} - ({p1_coins} + n2_1) = {total_value - p1_coins} - n2_1.")
    print("The game ends in a draw if V1 = V2, which means:")
    print(f"{p1_coins} + n2_1 = {total_value - p1_coins} - n2_1")
    print(f"2 * n2_1 = {total_value - 2 * p1_coins}")
    n2_1_threshold = (total_value - 2 * p1_coins)/2
    print(f"n2_1 = {n2_1_threshold}")
    print(f"\nSo, Player 1 wins if they collect more than {n2_1_threshold} of the {two_euro_coins} 2-euro coins.")
    print(f"Player 2 wins if Player 1 collects fewer than {n2_1_threshold} 2-euro coins.")
    print(f"It is a draw if Player 1 collects exactly {n2_1_threshold} 2-euro coins.\n")
    
    print("Step 3: Strategic Conclusion.")
    print("Because the total number of coins is odd, Player 1 gets the first and last pick.")
    print("This gives Player 1 more turns and more control over the game's outcome.")
    print("While a specific arrangement could lead to a draw, Player 2 can never force a win against an optimal Player 1.")
    print("Therefore, choosing to be the 1st player is always the better choice as it provides a non-losing strategy and the highest chance of winning.")

solve_coin_game()
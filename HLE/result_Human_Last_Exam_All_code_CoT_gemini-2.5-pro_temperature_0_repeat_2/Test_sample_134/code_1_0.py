def solve_coin_game():
    """
    Analyzes the coin game to determine if it's better to be the 1st or 2nd player.
    """
    # Step 1: Define the game parameters
    one_euro_coins = 136
    two_euro_coins = 87
    total_coins = one_euro_coins + two_euro_coins

    print(f"Total number of 1-euro coins: {one_euro_coins}")
    print(f"Total number of 2-euro coins: {two_euro_coins}")
    print(f"Total number of coins in the line: {total_coins}")
    print("-" * 30)

    # Step 2: Determine the number of coins per player
    # Since the total number of coins (223) is odd, the players take a different number of coins.
    # Player 1 (P1) makes moves 1, 3, 5, ..., 223.
    # Player 2 (P2) makes moves 2, 4, 6, ..., 222.
    p1_coins_count = (total_coins + 1) // 2
    p2_coins_count = total_coins // 2

    print(f"Player 1 picks {p1_coins_count} coins.")
    print(f"Player 2 picks {p2_coins_count} coins.")
    print("-" * 30)

    # Step 3: Formulate the scores and simplify the winning condition
    # Let n1_1 and n2_1 be the number of 1-euro and 2-euro coins P1 collects.
    # Let n1_2 and n2_2 be the number of 1-euro and 2-euro coins P2 collects.
    #
    # P1's score: V1 = 1 * n1_1 + 2 * n2_1
    # P2's score: V2 = 1 * n1_2 + 2 * n2_2
    #
    # We know:
    # n1_1 + n2_1 = p1_coins_count  => n1_1 = p1_coins_count - n2_1
    # n1_2 + n2_2 = p2_coins_count  => n1_2 = p2_coins_count - n2_2
    #
    # Substitute these into the score equations:
    # V1 = (p1_coins_count - n2_1) + 2 * n2_1 = p1_coins_count + n2_1
    # V2 = (p2_coins_count - n2_2) + 2 * n2_2 = p2_coins_count + n2_2
    
    print("Let n2_1 be the number of 2-euro coins for Player 1.")
    print("Let n2_2 be the number of 2-euro coins for Player 2.")
    print(f"Player 1's score can be expressed as: V1 = {p1_coins_count} + n2_1")
    print(f"Player 2's score can be expressed as: V2 = {p2_coins_count} + n2_2")
    print("-" * 30)

    # Step 4: Analyze the difference in scores
    # V1 - V2 = (p1_coins_count + n2_1) - (p2_coins_count + n2_2)
    # V1 - V2 = (p1_coins_count - p2_coins_count) + n2_1 - n2_2
    score_diff_constant = p1_coins_count - p2_coins_count
    
    print("The difference in scores is V1 - V2.")
    print(f"V1 - V2 = ({p1_coins_count} + n2_1) - ({p2_coins_count} + n2_2)")
    print(f"V1 - V2 = ({p1_coins_count} - {p2_coins_count}) + n2_1 - n2_2")
    print("The final equation for the score difference is:")
    print(f"V1 - V2 = {score_diff_constant} + n2_1 - n2_2")
    print("-" * 30)

    # Step 5: Interpret the result
    print("From the equation V1 - V2 = 1 + n2_1 - n2_2, we can determine the winning condition:")
    print("Player 1 wins if V1 > V2, which means 1 + n2_1 - n2_2 > 0, or n2_1 >= n2_2.")
    print("Players tie if V1 = V2, which means 1 + n2_1 - n2_2 = 0, or n2_2 = n2_1 + 1.")
    print("Player 2 wins if V2 > V1, which means 1 + n2_1 - n2_2 < 0, or n2_2 >= n2_1 + 2.")
    print("\nThe game is essentially a contest to collect more 2-euro coins.")
    print("\nPlayer 1 has a fundamental advantage because they get to pick one more coin than Player 2 (112 vs 111).")
    print("While some specific arrangements might allow Player 2 to win, for a random arrangement, the player who takes more coins is statistically favored.")
    print("Therefore, it is better to be the 1st player.")

solve_coin_game()
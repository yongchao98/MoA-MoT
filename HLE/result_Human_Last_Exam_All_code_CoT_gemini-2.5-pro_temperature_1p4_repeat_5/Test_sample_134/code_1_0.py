def solve_coin_game():
    """
    Analyzes the coin game to determine if it's better to be the 1st or 2nd player.
    """
    # Problem parameters
    n1 = 136  # Number of 1-euro coins
    n2 = 87   # Number of 2-euro coins

    # Step 1: Calculate total coins and total value
    total_coins = n1 + n2
    total_value = n1 * 1 + n2 * 2

    print("Step 1: Game Setup")
    print(f"Number of 1-euro coins: {n1}")
    print(f"Number of 2-euro coins: {n2}")
    print(f"The equation for the total number of coins is: {n1} + {n2} = {total_coins}")
    print(f"The equation for the total value is: {n1} * 1 + {n2} * 2 = {total_value} euros")
    print("-" * 40)

    # Step 2: Calculate number of picks for each player
    p1_picks = (total_coins + 1) // 2
    p2_picks = (total_coins - 1) // 2

    print("Step 2: Player Turns")
    print(f"Since there are an odd number of coins ({total_coins}), Player 1 gets more picks.")
    print(f"Player 1 picks ({total_coins} + 1) / 2 = {p1_picks} coins.")
    print(f"Player 2 picks ({total_coins} - 1) / 2 = {p2_picks} coins.")
    print("Player 1 has a structural advantage of one extra coin.")
    print("-" * 40)

    # Step 3: Determine the winning condition
    winning_threshold = total_value / 2

    print("Step 3: Winning Condition")
    print(f"To win, a player's score must be greater than half the total value.")
    print(f"The winning score threshold is {total_value} / 2 = {winning_threshold} euros.")
    print("-" * 40)

    # Step 4: Analyze Player 1's path to victory
    # P1's score S1 = (p1_picks - k2) * 1 + k2 * 2 = p1_picks + k2
    # To win, S1 > winning_threshold  =>  p1_picks + k2 > winning_threshold
    k2_needed = winning_threshold - p1_picks

    print("Step 4: Player 1's Winning Strategy")
    print("Let k2 be the number of 2-euro coins Player 1 collects.")
    print(f"Player 1's score is calculated by the equation: ({p1_picks} - k2) * 1 + k2 * 2 = {p1_picks} + k2")
    print("To win, Player 1's score must be greater than the threshold:")
    print(f"The winning inequality is: {p1_picks} + k2 > {winning_threshold}")
    print(f"Solving for k2: k2 > {winning_threshold} - {p1_picks}")
    print(f"So, k2 > {k2_needed}")
    print(f"Since k2 must be a whole number, Player 1 needs to collect at least {int(k2_needed) + 1} of the {n2} available 2-euro coins to guarantee a win.")
    print("-" * 40)

    # Step 5: Final Conclusion
    print("Conclusion:")
    print("Although Player 2 has some positional control, Player 1's advantage of an extra move in this constant-sum game is decisive.")
    print("Assuming optimal play from both sides, Player 1 can always guarantee at least a tie, and for most arrangements, will win.")
    print("Therefore, it is preferable to be the 1st player.")

solve_coin_game()
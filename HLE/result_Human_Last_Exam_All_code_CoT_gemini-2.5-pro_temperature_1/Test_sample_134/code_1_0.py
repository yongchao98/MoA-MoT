import math

def solve_coin_game():
    """
    Analyzes the coin game to determine if it's better to be the 1st or 2nd player.
    """
    n1 = 136  # Number of 1-euro coins
    n2 = 87   # Number of 2-euro coins
    
    n_total = n1 + n2
    total_value = n1 * 1 + n2 * 2
    
    # Step 1: Define the game structure based on odd number of coins
    n_black = math.ceil(n_total / 2)
    n_white = math.floor(n_total / 2)
    
    print("Step 1: Game Analysis")
    print(f"Total coins: {n_total} (an odd number)")
    print(f"Total value: {total_value} euros")
    print(f"Coloring positions reveals {n_black} 'black' spots and {n_white} 'white' spots.")
    print("This structure gives Player 2 control over which set of colored coins each player gets.")
    print("-" * 30)

    # Step 2: Define the winning condition
    print("Step 2: Winning Condition for Player 1")
    print("A detailed analysis shows that Player 1 wins if the following condition is met:")
    print("0 < S_B - S_W < 2 * c_max_ends")
    print("Where S_B is the sum of coins on black spots, S_W on white spots, and c_max_ends is the larger of the two end coins.")
    print("-" * 30)

    # Step 3: Calculate expected values to evaluate the condition on average
    
    # Expected value of a single random coin
    # E[coin] = 1 * P(1) + 2 * P(2) = 1 * (n1/n_total) + 2 * (n2/n_total)
    # E[coin] = (n1 + 2*n2) / n_total = total_value / n_total
    e_coin = total_value / n_total
    
    # Expected sum on black and white positions
    e_s_b = n_black * e_coin
    e_s_w = n_white * e_coin
    
    # Expected difference S_B - S_W
    e_diff_sw_sb = e_s_b - e_s_w
    
    # Calculate the expected value of the maximum of the two end coins
    # P(max=1) = P(c1=1 and c_end=1) = P(c1=1) * P(c_end=1 | c1=1)
    p_cmax_is_1 = (n1 / n_total) * ((n1 - 1) / (n_total - 1))
    p_cmax_is_2 = 1 - p_cmax_is_1
    e_cmax_ends = 1 * p_cmax_is_1 + 2 * p_cmax_is_2
    
    e_2_cmax_ends = 2 * e_cmax_ends
    
    print("Step 3: Calculating Expected Values")
    print(f"The expected sum on black positions E[S_B] is: {e_s_b:.2f}")
    print(f"The expected sum on white positions E[S_W] is: {e_s_w:.2f}")
    print(f"The expected difference E[S_B - S_W] is: {e_diff_sw_sb:.2f}")
    print(f"The expected value of the max of the two end coins E[c_max_ends] is: {e_cmax_ends:.2f}")
    print(f"The expected value of 2 * c_max_ends is: {e_2_cmax_ends:.2f}")
    print("-" * 30)

    # Step 4: Conclusion
    print("Step 4: Conclusion")
    print("We check if the average outcome satisfies Player 1's winning condition:")
    print(f"Is 0 < E[S_B - S_W] < E[2 * c_max_ends]?")
    print(f"Substituting the calculated values:")
    print(f"Is 0 < {e_diff_sw_sb:.2f} < {e_2_cmax_ends:.2f}?")
    
    if 0 < e_diff_sw_sb < e_2_cmax_ends:
        print("\nYes, the condition is met on average.")
        print("This means that in a typical random arrangement, Player 1 has a winning strategy.")
        print("Therefore, it is preferable to be the 1st player.")
    else:
        print("\nNo, the condition is not met on average, suggesting Player 2 has the advantage.")

solve_coin_game()
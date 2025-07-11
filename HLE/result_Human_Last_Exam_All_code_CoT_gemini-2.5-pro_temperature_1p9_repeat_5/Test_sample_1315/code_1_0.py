def solve_tichu_problem():
    """
    Calculates the maximal possible value of X-Y in a single Tichu round
    under the given conditions.
    """

    # --- Step 1: Define point values in Tichu ---
    total_card_points = 100
    grand_tichu_success = 200
    grand_tichu_failure = -200

    # --- Step 2: Determine the optimal scenario for each component ---

    # To maximize X-Y, the winning team (Team A) must maximize their score
    # and minimize the losing team's (Team B) score.

    # Maximize Team A's card points. This happens if they capture all point cards.
    # This is possible if the last player to go out is from Team B, and they
    # give their captured points to the first player (from Team A).
    card_points_A = total_card_points
    card_points_B = total_card_points - card_points_A

    # Maximize Team A's Tichu points and minimize Team B's.
    # Team A: One player makes a successful Grand Tichu call.
    # Team B: Both players make a failed Grand Tichu call.
    tichu_points_A = grand_tichu_success
    num_failed_tichus_B = 2
    tichu_points_B = num_failed_tichus_B * grand_tichu_failure

    # --- Step 3: Calculate the scores X and Y ---
    X = card_points_A + tichu_points_A
    Y = card_points_B + tichu_points_B

    # --- Step 4: Calculate the final difference and print the results ---
    difference = X - Y

    print("To achieve the maximal difference (X-Y), the following must happen:")
    print("1. The winning team (Team A) captures all 100 card points.")
    print("2. One player from Team A makes a successful 'Grand Tichu' call.")
    print("3. Both players from the losing team (Team B) make failed 'Grand Tichu' calls.")
    print("\nCalculating the scores:")
    
    print(f"Winning team's score (X) = {card_points_A} (card points) + {tichu_points_A} (Grand Tichu) = {X}")
    print(f"Losing team's score (Y) = {card_points_B} (card points) + {tichu_points_B} (2 failed Grand Tichus) = {Y}")
    
    print("\nCalculating the maximal difference (X - Y):")
    # To satisfy the prompt's request to show each number in the final equation.
    if Y < 0:
        print(f"X - Y = {X} - ({Y}) = {difference}")
    else:
        print(f"X - Y = {X} - {Y} = {difference}")
    
    # Final answer in the required format
    print(f"\n<<<{difference}>>>")

solve_tichu_problem()
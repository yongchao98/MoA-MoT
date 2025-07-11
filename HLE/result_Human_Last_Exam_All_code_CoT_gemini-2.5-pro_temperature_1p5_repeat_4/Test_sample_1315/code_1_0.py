def solve_tichu_scenario():
    """
    Calculates the maximal possible score difference (X-Y) in a Tichu round
    where the winning team does not finish 1st and 2nd.

    The scenario assumes:
    - Team A (winning team) players finish 1st and 3rd.
    - Team B (losing team) players finish 2nd and last.
    """

    # 1. Points for going out first
    # Team A's player finishes first, getting 100 points.
    first_place_bonus = 100
    print(f"Team A gets a bonus of {first_place_bonus} for finishing first.")

    # 2. Points from Tichu calls
    # To maximize the difference, we assume the most extreme valid Tichu calls.
    # Team A's first player makes a successful Grand Tichu.
    tichu_A_score = 200
    print(f"Team A scores {tichu_A_score} from a successful Grand Tichu call.")
    # Both players on Team B attempt a Tichu call and fail.
    tichu_B_score = -100 * 2
    print(f"Team B scores {tichu_B_score} from two failed Tichu calls.")

    # 3. Points from captured cards
    # The total value of point cards in the deck is 100.
    # The player from Team B who finished last gives their captured tricks
    # and their remaining hand cards to Team A.
    # To maximize the difference, we assume the other player from Team B captured 0 points.
    # Therefore, all 100 card points go to Team A.
    card_points_A = 100
    card_points_B = 0
    print(f"Team A collects all {card_points_A} available card points.")
    print(f"Team B collects {card_points_B} card points.")

    # 4. Calculate final scores X and Y
    X = first_place_bonus + tichu_A_score + card_points_A
    Y = tichu_B_score + card_points_B
    
    print("\n--- Final Score Calculation ---")
    print(f"Winning team's score (X) = {first_place_bonus} (1st place) + {card_points_A} (cards) + {tichu_A_score} (Grand Tichu) = {X}")
    print(f"Losing team's score (Y) = {card_points_B} (cards) + ({tichu_B_score}) (failed Tichus) = {Y}")

    # 5. Calculate the difference X-Y
    difference = X - Y
    print("\n--- Maximal Difference (X-Y) ---")
    print(f"The final equation is: {X} - ({Y}) = {difference}")

solve_tichu_scenario()
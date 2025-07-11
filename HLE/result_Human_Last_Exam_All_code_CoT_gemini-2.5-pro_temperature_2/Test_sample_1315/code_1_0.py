def solve_tichu_problem():
    """
    Calculates the maximal possible score difference (X-Y) in a round of Tichu
    under the specified conditions.

    The scenario for maximal difference is as follows:
    1. Card Points: The winning team (Team A) gets all 100 card points. The
       losing team (Team B) gets 0.
    2. Calls:
       - One player on Team A makes a successful "Grand Tichu" call (+200 points).
       - Both players on Team B make unsuccessful "Grand Tichu" calls (-200 points each).
    3. Condition Check: This requires a player from Team A to go out first. Since
       the winning team does not go out first and second, the order could be
       A, B, A, B, which is a valid round finish.
    """

    # --- Points for the Winning Team (X) ---

    # Card points: Winning team collects all 100 points.
    winning_team_card_points = 100

    # Call points: One player makes a successful Grand Tichu call.
    winning_team_call_points = 200

    # Total score for the winning team (X).
    X = winning_team_card_points + winning_team_call_points

    # --- Points for the Losing Team (Y) ---

    # Card points: Losing team collects 0 points from cards.
    losing_team_card_points = 0

    # Call points: Both players make a Grand Tichu call and fail.
    losing_team_failed_call_1 = -200
    losing_team_failed_call_2 = -200
    losing_team_call_points = losing_team_failed_call_1 + losing_team_failed_call_2

    # Total score for the losing team (Y).
    Y = losing_team_card_points + losing_team_call_points

    # --- Calculate the Maximal Difference (X - Y) ---
    max_difference = X - Y

    # Print the breakdown of the calculation as requested.
    print(f"To achieve the maximal difference (X-Y):")
    print(f"The winning team's score (X) must be maximized.")
    print(f"X = (Card Points) + (Call Points) = {winning_team_card_points} + {winning_team_call_points} = {X}")
    print("\n")
    print(f"The losing team's score (Y) must be minimized.")
    print(f"Y = (Card Points) + (Call Points) = {losing_team_card_points} + ({losing_team_failed_call_1} + {losing_team_failed_call_2}) = {Y}")
    print("\n")
    print("The final calculation for the maximal difference is:")
    print(f"{X} - ({Y}) = {max_difference}")

solve_tichu_problem()
<<<700>>>
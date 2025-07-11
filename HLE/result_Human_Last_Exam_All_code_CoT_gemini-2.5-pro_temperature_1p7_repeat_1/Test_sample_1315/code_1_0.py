def solve_tichu_problem():
    """
    Calculates the maximal possible score difference (X-Y) in a Tichu round
    under the specified constraints.

    - X: Winning team's score
    - Y: Losing team's score
    - Constraint: Winning team does not finish 1st and 2nd.
    """

    # --- Step 1: Scenario Definition ---
    # The optimal scenario for maximizing X-Y is as follows:
    # Teams: Team A (Winner) vs. Team B (Loser)
    # Finishing Order: A1, B1, A2, B2
    # This order satisfies the constraint and allows for the last player's
    # points (B2) to be transferred to the first player (A1).

    # Points from Tichu calls
    GRAND_TICHU_SUCCESS = 200
    GRAND_TICHU_FAILURE = -200

    # Total points from cards in the deck
    TOTAL_CARD_POINTS = 100

    # --- Step 2: Calculate Card Point Difference ---
    # To maximize the winner's score, we assume they collect all card points.
    # This is possible if player B1 collects 0 points, and player B2's points
    # are transferred to Team A.
    card_points_winner = TOTAL_CARD_POINTS
    card_points_loser = 0

    # --- Step 3: Calculate Tichu Point Difference ---
    # Winner's Tichu Points: Player A1 succeeds in a "Grand Tichu".
    tichu_points_winner = GRAND_TICHU_SUCCESS

    # Loser's Tichu Points: Players B1 and B2 both fail their "Grand Tichu" calls.
    tichu_points_loser = GRAND_TICHU_FAILURE + GRAND_TICHU_FAILURE

    # --- Step 4: Calculate Final Scores (X and Y) ---
    # X is the winning team's total score.
    X = card_points_winner + tichu_points_winner

    # Y is the losing team's total score.
    Y = card_points_loser + tichu_points_loser

    # --- Step 5: Calculate and Display the Final Difference ---
    difference = X - Y

    print("Calculation for the maximal difference X-Y in a Tichu round:")
    print("-" * 60)
    print(f"Winning Team Score (X):")
    print(f"  Card Points: {card_points_winner}")
    print(f"  Tichu Points: {tichu_points_winner} (Successful Grand Tichu)")
    print(f"  X = {card_points_winner} + {tichu_points_winner} = {X}")
    print("-" * 60)
    print(f"Losing Team Score (Y):")
    print(f"  Card Points: {card_points_loser}")
    print(f"  Tichu Points: {tichu_points_loser} (Two failed Grand Tichus)")
    print(f"  Y = {card_points_loser} + {tichu_points_loser} = {Y}")
    print("-" * 60)
    print(f"The final maximal difference (X - Y) is:")
    print(f"{X} - ({Y}) = {difference}")


solve_tichu_problem()
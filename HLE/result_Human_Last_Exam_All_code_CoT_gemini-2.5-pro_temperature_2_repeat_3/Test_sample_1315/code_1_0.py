def solve_tichu_problem():
    """
    Calculates the maximal possible value of X-Y in a Tichu game round
    where the winning team (score X) does not go out first and second.
    """

    # --- Game Scoring Parameters ---
    total_card_points = 100
    grand_tichu_bonus = 200
    grand_tichu_penalty = -200

    print("To find the maximal value of X-Y, we construct the most extreme scoring scenario possible.")
    print("The winning team is Team A (score X), and the losing team is Team B (score Y).")
    print("The finish order must be 1st/3rd vs 2nd/4th.\n")

    # --- Maximizing X (Score of Winning Team A) ---
    print("--- Calculating the Winning Team's Score (X) ---")
    # Team A captures all available points from cards
    team_A_card_points = total_card_points
    print(f"To maximize their score, Team A captures all card points: {team_A_card_points}")

    # A player on Team A successfully calls "Grand Tichu" by going out first
    team_A_bonus = grand_tichu_bonus
    print(f"A player on Team A makes a successful Grand Tichu call: +{team_A_bonus} points")

    X = team_A_card_points + team_A_bonus
    print(f"X = {team_A_card_points} (cards) + {team_A_bonus} (call) = {X}\n")

    # --- Minimizing Y (Score of Losing Team B) ---
    print("--- Calculating the Losing Team's Score (Y) ---")
    # Team B gets no points from cards
    team_B_card_points = 0
    print(f"Team B consequently captures no card points: {team_B_card_points}")

    # Both players on Team B call "Grand Tichu" and fail
    team_B_penalty_1 = grand_tichu_penalty
    team_B_penalty_2 = grand_tichu_penalty
    total_team_B_penalty = team_B_penalty_1 + team_B_penalty_2
    print(f"To minimize their score, both players on Team B make failed Grand Tichu calls:")
    print(f"  Player 1 penalty: {team_B_penalty_1} points")
    print(f"  Player 2 penalty: {team_B_penalty_2} points")


    Y = team_B_card_points + total_team_B_penalty
    print(f"Y = {team_B_card_points} (cards) + ({team_B_penalty_1}) (call 1) + ({team_B_penalty_2}) (call 2) = {Y}\n")

    # --- Calculating the Final Difference (X-Y) ---
    print("--- Calculating the Final Difference (X - Y) ---")
    difference = X - Y
    print(f"The maximal difference is X - Y.")
    print(f"X - Y = {X} - ({Y}) = {difference}")


solve_tichu_problem()
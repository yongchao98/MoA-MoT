def calculate_max_tichu_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a single Tichu round
    under the constraint that the winning team does not finish 1st and 2nd.
    """

    # Step 1: Maximize the card point difference.
    # To maximize the difference, the winning team (X) gets all 100 card points.
    # The losing team (Y) gets 0 points.
    card_points_X = 100
    card_points_Y = 0
    print("--- Card Point Calculation ---")
    print(f"To maximize the difference, we assume the winning team X gets all card points.")
    print(f"Card points for winning team X: {card_points_X}")
    print(f"Card points for losing team Y: {card_points_Y}\n")

    # Step 2: Maximize the Tichu call score difference.
    # The winning team (X) gets the max possible bonus.
    # The losing team (Y) gets the max possible penalty.
    
    # Team X: A player makes a successful "Grand Tichu" call (+200 points).
    tichu_score_X = 200
    
    # Team Y: Both players make "Grand Tichu" calls and fail (-200 points each).
    tichu_score_Y = -200 + -200

    print("--- Tichu Call Score Calculation ---")
    print("For winning team X:")
    print(f"  A successful 'Grand Tichu' call results in a bonus of: +{tichu_score_X}")
    print("For losing team Y:")
    print(f"  Two failed 'Grand Tichu' calls result in a penalty of: -200 and -200, for a total of {tichu_score_Y}\n")

    # Step 3: Calculate the total scores for each team.
    # Total Score = Card Points + Tichu Score
    X = card_points_X + tichu_score_X
    Y = card_points_Y + tichu_score_Y
    
    print("--- Final Score Calculation ---")
    print("Score for winning team X = (Card Points) + (Tichu Score)")
    print(f"X = {card_points_X} + {tichu_score_X} = {X}")
    print("Score for losing team Y = (Card Points) + (Tichu Score)")
    print(f"Y = {card_points_Y} + ({tichu_score_Y}) = {Y}\n")

    # Step 4: Calculate the final difference X - Y.
    difference = X - Y
    
    print("--- Maximal Difference Calculation ---")
    print("The maximal possible value of X - Y is calculated as:")
    print(f"{X} - ({Y}) = {difference}")

if __name__ == '__main__':
    calculate_max_tichu_difference()
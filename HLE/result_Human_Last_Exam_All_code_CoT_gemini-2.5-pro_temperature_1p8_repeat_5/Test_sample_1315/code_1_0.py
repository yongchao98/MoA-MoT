def solve_tichu_problem():
    """
    Calculates the maximal possible score difference (X-Y) in a single Tichu round
    where the winning team (X) does not go out first and second.
    """
    # Step 1: Explain the methodology
    print("To find the maximal possible value of X-Y, we must construct a scenario that maximizes the winning team's score (X) and minimizes the losing team's score (Y).")
    print("A team's score is the sum of their captured card points and Tichu call bonuses.\n")
    print("The plan is as follows:")
    print("1. Maximize the winning team's card points.")
    print("2. Maximize the Tichu bonus difference between the teams.")
    print("3. Calculate the final scores and their difference.\n")

    # Step 2: Define the optimal card point distribution
    points_kings = 40
    points_tens = 40
    points_fives = 20
    points_dragon = 25
    points_phoenix = -25

    # The winning team (X) captures all positive point cards and the Dragon.
    card_points_X = points_kings + points_tens + points_fives + points_dragon
    # The losing team (Y) captures the Phoenix.
    card_points_Y = points_phoenix

    print("--- Maximizing Card Points ---")
    print(f"The winning team (X) captures all Kings ({points_kings}), Tens ({points_tens}), Fives ({points_fives}), and the Dragon ({points_dragon}).")
    print(f"Card Points for Team X = {points_kings} + {points_tens} + {points_fives} + {points_dragon} = {card_points_X}")
    print(f"The losing team (Y) gets the trick with the Phoenix.")
    print(f"Card Points for Team Y = {points_phoenix}\n")

    # Step 3: Define the optimal Tichu call outcomes
    # Winning team succeeds in a Grand Tichu
    tichu_bonus_X = 200
    # Losing team fails a Grand Tichu
    tichu_bonus_Y = -200

    print("--- Maximizing Tichu Bonus Difference ---")
    print("A player on Team X successfully calls a 'Grand Tichu', earning +200 points.")
    print(f"Tichu Bonus for Team X = {tichu_bonus_X}")
    print("A player on Team Y also calls 'Grand Tichu' but fails (as Team X went out first), incurring a -200 point penalty.")
    print(f"Tichu Bonus for Team Y = {tichu_bonus_Y}\n")

    # Step 4: Calculate the total scores for each team
    X = card_points_X + tichu_bonus_X
    Y = card_points_Y + tichu_bonus_Y

    print("--- Final Score Calculation ---")
    print(f"Total Score for Winning Team (X) = Card Points + Tichu Bonus = {card_points_X} + {tichu_bonus_X} = {X}")
    print(f"Total Score for Losing Team (Y) = Card Points + Tichu Bonus = {card_points_Y} + {tichu_bonus_Y} = {Y}\n")

    # Step 5: Calculate and print the final difference
    difference = X - Y
    
    print("--- Maximal Difference (X - Y) ---")
    print("The final equation for the difference is:")
    # Using f-string to ensure all numbers in the final equation are shown
    print(f"({card_points_X} + {tichu_bonus_X}) - ({card_points_Y} + ({tichu_bonus_Y})) = {X} - ({Y}) = {difference}")

if __name__ == "__main__":
    solve_tichu_problem()
def solve_tichu_max_difference():
    """
    Calculates the maximal score difference in a game of Tichu under specific conditions.
    """
    # Step 1: Define the optimal point distribution for the scenario.

    # To maximize X, the winning team gets all card points and a Grand Tichu bonus.
    # The last player to go out is on the losing team. They give their captured tricks
    # to the team that went out first (the winners). Thus, the winners can get all 100 card points.
    card_points_winner = 100
    # The first player out (a winner) made a successful Grand Tichu.
    tichu_bonus_winner = 200

    # To minimize Y, the losing team gets no card points and maximum penalties.
    card_points_loser = 0
    # Both players on the losing team unsuccessfully call Grand Tichu (-200 each).
    tichu_penalty_loser = -400

    # Step 2: Calculate the scores X (winner) and Y (loser).
    X = card_points_winner + tichu_bonus_winner
    Y = card_points_loser + tichu_penalty_loser

    # Step 3: Calculate the difference and print the detailed equation.
    difference = X - Y

    print("This script calculates the maximal score difference (X - Y) in a round of Tichu.")
    print("The winning team (score X) is assumed not to have a 'one-two' victory.\n")

    print("--- Calculating the Scores ---")
    print(f"Winning Team's Score (X) = (Card Points) + (Tichu Bonus)")
    print(f"X = {card_points_winner} + {tichu_bonus_winner} = {X}\n")

    print(f"Losing Team's Score (Y) = (Card Points) + (Tichu Penalty)")
    print(f"Y = {card_points_loser} + ({tichu_penalty_loser}) = {Y}\n")

    print("--- Calculating the Final Difference (X - Y) ---")
    print("The final equation is built from the components of each team's score:")
    print(f"X - Y = ({card_points_winner} + {tichu_bonus_winner}) - ({card_points_loser} + ({tichu_penalty_loser}))")
    print(f"X - Y = {X} - ({Y})")
    print(f"X - Y = {difference}")


solve_tichu_max_difference()

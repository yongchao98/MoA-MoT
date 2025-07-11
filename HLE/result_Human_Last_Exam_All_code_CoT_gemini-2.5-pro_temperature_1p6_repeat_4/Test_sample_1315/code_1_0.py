import sys

def solve_tichu_max_score_difference():
    """
    Calculates the maximal possible value of X-Y in a single Tichu round,
    where X is the winning team's score and Y is the losing team's score,
    given the winning team does not go out first and second.
    """

    # --- Setup of Game Values ---

    # Total points from cards (Kings, Tens, Fives) in a round.
    total_card_points = 100
    # Point value for a successful or failed Grand Tichu call.
    grand_tichu_value = 200

    # --- Step 1: Maximize the Card Point Difference ---

    # To maximize the difference, the winning team (X) must acquire all card points.
    # This is possible because the player finishing last gives their tricks and hand to the winners.
    winning_team_max_card_points = 100
    losing_team_min_card_points = 0

    # --- Step 2: Maximize the Tichu Call Point Difference ---

    # The winning team's score is maximized with a successful Grand Tichu call by the first player out.
    winning_team_tichu_points = grand_tichu_value

    # The losing team's score is minimized when both players fail a Grand Tichu call.
    losing_player1_failed_tichu = -grand_tichu_value
    losing_player2_failed_tichu = -grand_tichu_value
    losing_team_tichu_points = losing_player1_failed_tichu + losing_player2_failed_tichu

    # --- Step 3: Calculate the Total Scores (X and Y) ---

    # X = Winning team's total score
    X = winning_team_max_card_points + winning_team_tichu_points

    # Y = Losing team's total score
    Y = losing_team_min_card_points + losing_team_tichu_points

    # --- Step 4: Calculate the Final Difference and Print the Derivation ---

    print("This script calculates the maximal score difference (X-Y) in a Tichu round without a 1-2 finish.\n")
    print("Let X be the winning team's score and Y be the losing team's score.\n")

    print("--- SCORE CALCULATION ---")
    print("1. Maximizing Card Points:")
    print(f"The winning team acquires all {winning_team_max_card_points} card points. The losing team gets {losing_team_min_card_points}.\n")

    print("2. Maximizing Tichu Call Points:")
    print(f"The winning team makes a successful Grand Tichu for +{winning_team_tichu_points} points.")
    print(f"Both players on the losing team fail Grand Tichu calls, for a total of {losing_team_tichu_points} points.\n")

    print("--- FINAL SCORES ---")
    print(f"Winning Team Score (X) = {winning_team_max_card_points} (cards) + {winning_team_tichu_points} (Tichu) = {X}")
    print(f"Losing Team Score (Y) = {losing_team_min_card_points} (cards) + ({losing_team_tichu_points}) (Tichu) = {Y}\n")

    final_difference = X - Y
    print("--- MAXIMAL DIFFERENCE (X - Y) ---")
    # The final equation is printed here as requested.
    print(f"The final calculation is: {X} - ({Y}) = {final_difference}")

    # Output the final answer in the specified format
    # Redirecting to stderr to keep the final answer clean for parsing
    print(f"\n<<<{final_difference}>>>", file=sys.stderr)

solve_tichu_max_score_difference()
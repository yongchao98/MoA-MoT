def solve_tichu_problem():
    """
    Calculates the maximal possible value of X-Y in a single Tichu round
    where the winning team (score X) does not go out first and second.
    """

    # Step 1: Maximize the Tichu/Grand Tichu bonus difference.
    # To maximize the difference, the winning team (A) gets the highest possible bonus,
    # and the losing team (B) gets the biggest penalty.
    # Team A: Successful "Grand Tichu" call -> +200 points.
    # Team B: Failed "Tichu" call -> -100 points.
    tichu_bonus_A = 200
    tichu_bonus_B = -100

    # Step 2: Maximize the card point difference.
    # The total card points in the deck is 100.
    # To maximize the difference, Team A should get all the points.
    # This is possible if the last player is on Team B, and the other Team B player
    # takes no tricks with point cards.
    # In this scenario, all points taken by Team A players, plus all points taken by
    # the last player from Team B, and all points in the last player's hand go to Team A.
    card_points_A = 100
    card_points_B = 0

    # Step 3: Calculate the final scores for each team.
    # X = Winning Team Score, Y = Losing Team Score
    X = card_points_A + tichu_bonus_A
    Y = card_points_B + tichu_bonus_B

    # Step 4: Calculate the final difference.
    difference = X - Y

    # Output the logic and the final calculation breakdown.
    print("To find the maximal difference X-Y, we need to maximize the winning team's score (X) and minimize the losing team's score (Y).")
    print("\nThis is achieved under the following optimal (but possible) scenario:")
    print("1. A player on the winning team (Team A) successfully calls 'Grand Tichu', earning +200 points.")
    print("2. A player on the losing team (Team B) fails a 'Tichu' call, incurring a -100 point penalty.")
    print("3. The finishing order is A, B, A, with the last player being on Team B. This arrangement ensures that all 100 card points are awarded to Team A.")

    print("\n--- Final Calculation ---")
    print(f"Winning team's score (X) = (Card Points) + (Grand Tichu Bonus)")
    print(f"X = {card_points_A} + {tichu_bonus_A} = {X}")
    print(f"\nLosing team's score (Y) = (Card Points) + (Tichu Penalty)")
    print(f"Y = {card_points_B} + ({tichu_bonus_B}) = {Y}")
    print("\n--- Maximal Difference ---")
    print(f"X - Y = {X} - ({Y}) = {difference}")

solve_tichu_problem()
<<<400>>>
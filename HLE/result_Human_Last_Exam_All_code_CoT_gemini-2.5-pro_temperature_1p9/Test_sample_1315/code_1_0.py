def solve_tichu_max_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a single Tichu round,
    given the winning team does not finish first and second.
    """

    # 1. Card Points Analysis
    # Total points from cards in the deck (5s, 10s, Kings, Dragon, Phoenix).
    total_card_points = 100

    # To maximize (C_win - C_lose) where C_win + C_lose = 100, we must maximize C_win.
    # This is achieved if the winning team collects all 100 card points.
    # A player from the winning team finishes first, and a player from the losing team
    # finishes last. The last player's tricks are given to the first.
    # If the other losing team player takes 0 points, the winning team gets all 100 points.
    c_win = 100
    c_lose = 0
    
    print("--- Score Calculation Breakdown ---")
    print("\n1. Card Point Distribution:")
    print(f"To maximize the difference, the winning team must get all card points.")
    print(f"Winning Team Card Points (C_win): {c_win}")
    print(f"Losing Team Card Points (C_lose): {c_lose}")

    # 2. Tichu Call Analysis
    # To maximize (T_win - T_lose), we must maximize T_win and minimize T_lose.
    # T_win is maximized if a player on the winning team succeeds at a "Grand Tichu" (+200).
    # This is possible because a player from their team finishes first.
    t_win = 200

    # T_lose is minimized if both players on the losing team fail "Grand Tichu" calls (-200 each).
    # This is guaranteed as they do not finish first.
    t_lose = -200 + -200

    print("\n2. Tichu Call Point Distribution:")
    print("Winning Team: One player calls and wins a 'Grand Tichu'.")
    print(f"Winning Team Call Points (T_win): {t_win}")
    print("Losing Team: Both players call 'Grand Tichu' and fail.")
    print(f"Losing Team Call Points (T_lose): {-200} + {-200} = {t_lose}")

    # 3. Final Score Calculation
    # X = Winning team's total score
    # Y = Losing team's total score
    X = c_win + t_win
    Y = c_lose + t_lose
    
    difference = X - Y
    
    print("\n--- Final Score ---")
    print(f"Winning Team Score (X) = {c_win} + {t_win} = {X}")
    print(f"Losing Team Score (Y) = {c_lose} + ({t_lose}) = {Y}")

    print("\n--- Maximal Difference (X - Y) ---")
    print(f"The equation is: {X} - ({Y})")
    print(f"Maximal Possible Value of X - Y = {difference}")


solve_tichu_max_difference()
<<<700>>>
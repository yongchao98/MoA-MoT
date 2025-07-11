def calculate_max_tichu_score_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a single Tichu round
    where the winning team does not go out first and second.
    """

    # 1. Calculate the maximum difference from card points.
    # The total points from cards in a deck is 100.
    # Points are: Kings (10), Tens (10), Fives (5), Dragon (25), Phoenix (-25).
    # Total positive points = 4*10 (Kings) + 4*10 (Tens) + 4*5 (Fives) + 25 (Dragon) = 125.
    # Total negative points = -25 (Phoenix).
    # To maximize the winning team's score (X_cards) and minimize the losing team's (Y_cards),
    # the winning team must get all the positive cards and the losing team the negative one.
    # This is possible in a regular round (not a 1-2 finish).
    # For example, if the last player is on the losing team, they give their captured
    # tricks to the winning team. If the other player on the losing team has the Phoenix
    # and nothing else, the score split becomes 125 and -25.
    max_X_cards = 125
    min_Y_cards = -25

    # 2. Calculate the maximum difference from bonus/penalty points.
    # To maximize the winning team's score (X), their player who goes out first
    # must have made a successful "Grand Tichu" call.
    max_X_bonus = 200  # Successful Grand Tichu call for the winning team

    # To minimize the losing team's score (Y), both players on that team must
    # have made a "Tichu" call and failed (because the other team's player went out first).
    min_Y_bonus = -100 + -100  # Two failed Tichu calls for the losing team

    # 3. Calculate the total scores X and Y.
    max_X = max_X_cards + max_X_bonus
    min_Y = min_Y_cards + min_Y_bonus

    # 4. Calculate the final difference X - Y.
    result = max_X - min_Y

    # Print the final equation with all components
    print("Maximum score for winning team (X):")
    print(f"  Card Points: {max_X_cards}")
    print(f"  Bonus Points: {max_X_bonus} (Successful Grand Tichu)")
    print(f"  Total X = {max_X}\n")

    print("Minimum score for losing team (Y):")
    print(f"  Card Points: {min_Y_cards}")
    print(f"  Penalty Points: {min_Y_bonus} (Two failed Tichus)")
    print(f"  Total Y = {min_Y}\n")

    print("Calculation of the maximal difference (X - Y):")
    # The format X - Y needs careful formatting for the negative Y value.
    print(f"({max_X_cards} + {max_X_bonus}) - ({min_Y_cards} + {min_Y_bonus}) = {max_X} - ({min_Y}) = {result}")

calculate_max_tichu_score_difference()
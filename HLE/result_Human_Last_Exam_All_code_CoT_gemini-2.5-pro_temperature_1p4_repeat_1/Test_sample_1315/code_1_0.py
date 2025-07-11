def solve_tichu_max_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a standard
    Tichu round, given the winning team does not finish first and second.
    """

    # 1. Deconstruct Card Scores
    # The total points available from cards in the deck is 100.
    # To maximize the difference, we assume the winning team (X) captures all
    # possible points, and the losing team (Y) captures none.
    # This is theoretically possible through the game's card transfer rules.
    x_cards = 100
    y_cards = 0

    # 2. Deconstruct Tichu Bonuses
    # "Grand Tichu" calls provide the largest score swing (+200 for success, -200 for failure).
    # We construct the optimal scenario for maximizing the bonus difference.

    # Winning team's bonus (X_bonus):
    # One player on the winning team declares "Grand Tichu" and goes out first.
    # This succeeds, earning the team 200 points.
    x_bonus = 200

    # Losing team's bonus (Y_bonus):
    # Both players on the losing team declare "Grand Tichu". Since neither can
    # finish first, both calls fail, resulting in a 200-point penalty for each.
    y_bonus = -200 - 200

    # 3. Calculate Final Scores (X and Y)
    # The total score for each team is the sum of their card points and bonuses.
    X = x_cards + x_bonus
    Y = y_cards + y_bonus

    # 4. Calculate the Maximal Difference
    difference = X - Y

    # Print the breakdown of the calculation as requested
    print("The equation for the final score of the winning team (X) is:")
    print(f"X = (Card Points) + (Tichu Bonus)")
    print(f"{x_cards} + {x_bonus} = {X}")
    print("\nThe equation for the final score of the losing team (Y) is:")
    print(f"Y = (Card Points) + (Tichu Bonus)")
    print(f"{y_cards} + ({y_bonus}) = {Y}")
    print("\nThe final equation for the maximal difference (X - Y) is:")
    print(f"({X}) - ({Y}) = {difference}")

solve_tichu_max_difference()

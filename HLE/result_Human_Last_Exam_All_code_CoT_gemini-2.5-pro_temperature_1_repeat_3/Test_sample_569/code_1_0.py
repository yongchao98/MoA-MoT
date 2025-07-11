def calculate_pinochle_score():
    """
    Calculates the total points for a Pinochle hand based on the provided scenario.
    """
    # 1. Calculate the declarer's (my) meld points
    my_run_points = 150
    # 8 aces = 2 sets of "aces around" (100 points per set)
    my_aces_points = 2 * 100
    my_total_meld = my_run_points + my_aces_points

    # 2. Calculate the partner's meld points
    partner_pinochle_points = 40
    # Both 9s of trump = 2 * 10 points
    partner_nines_points = 2 * 10
    partner_total_meld = partner_pinochle_points + partner_nines_points

    # 3. Calculate points from tricks
    # With a run, 8 aces, and "perfect play", all tricks will be won.
    # Total points from counter cards in the deck is 240.
    counter_points = 240
    # A 10-point bonus is awarded for winning the last trick.
    last_trick_bonus = 10
    total_trick_points = counter_points + last_trick_bonus

    # 4. Sum all points for the final score
    final_score = my_total_meld + partner_total_meld + total_trick_points

    # Print the detailed breakdown of the final score equation
    print("The total score is calculated by adding all meld and trick points:")
    print(f"({my_run_points} + {my_aces_points}) [My Meld] + "
          f"({partner_pinochle_points} + {partner_nines_points}) [Partner's Meld] + "
          f"({counter_points} + {last_trick_bonus}) [Trick Points] = {final_score}")

calculate_pinochle_score()
<<<660>>>
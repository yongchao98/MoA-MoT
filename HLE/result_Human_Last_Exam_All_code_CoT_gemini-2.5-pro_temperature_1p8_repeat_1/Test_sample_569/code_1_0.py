def calculate_pinochle_score():
    """
    Calculates the total score for a Pinochle hand based on a specific scenario.
    """
    # 1. My Meld
    my_run = 15
    my_aces = 100
    my_total_meld = my_run + my_aces
    print(f"My meld points: {my_run} (run) + {my_aces} (8 aces) = {my_total_meld} points.")

    # 2. Partner's Meld
    partner_pinochle = 4
    partner_nines = 10 * 2  # Two 9s of trump (dix)
    partner_total_meld = partner_pinochle + partner_nines
    print(f"Partner's meld points: {partner_pinochle} (pinochle) + {partner_nines} (two 9s) = {partner_total_meld} points.")

    # 3. Total Team Meld
    total_meld = my_total_meld + partner_total_meld
    print(f"Total team meld: {my_total_meld} + {partner_total_meld} = {total_meld} points.")
    print("-" * 30)

    # 4. Points from Play
    # Point values: A=11, 10=10, K=4, Q=3, J=2, 9=0. There are 8 of each.
    points_from_cards = (8 * 11) + (8 * 10) + (8 * 4) + (8 * 3) + (8 * 2) + (8 * 0)
    last_trick_bonus = 10
    
    # "Playing perfectly" with this hand means winning all tricks.
    total_play_points = points_from_cards + last_trick_bonus
    print("With a run and 8 aces, playing perfectly means winning all tricks.")
    print(f"Points from play: {points_from_cards} (all cards) + {last_trick_bonus} (last trick bonus) = {total_play_points} points.")
    print("-" * 30)

    # 5. Final Score
    final_score = total_meld + total_play_points
    print(f"Total score for the hand: {total_meld} (meld) + {total_play_points} (play) = {final_score} points.")

calculate_pinochle_score()
<<<389>>>
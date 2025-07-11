def calculate_pinochle_score():
    """
    Calculates the total points for a specific hand of Pinochle based on the user's description.
    """

    # --- Meld Points ---
    # My meld
    my_run_points = 15  # A, 10, K, Q, J of trump
    my_aces_around_points = 10  # One of each Ace
    my_num_aces_around = 2 # Having 8 aces provides two sets
    my_total_aces_points = my_aces_around_points * my_num_aces_around
    my_total_meld = my_run_points + my_total_aces_points

    # Partner's meld
    partner_pinochle_points = 4 # Q of Spades and J of Diamonds
    partner_dix_points = 10 # 9 of trump
    partner_num_dix = 2
    partner_total_dix_points = partner_dix_points * partner_num_dix
    partner_total_meld = partner_pinochle_points + partner_total_dix_points
    
    # --- Points from Play (Tricks) ---
    # In a 48-card deck, Aces, 10s, and Kings are counters (1 point each).
    # There are 8 of each (2x per suit). Total = 8+8+8 = 24.
    counter_card_points = 24
    last_trick_bonus = 1
    # Given the described hand, playing perfectly guarantees winning all tricks.
    # Therefore, the team collects all counter points and the last trick bonus.
    total_play_points = counter_card_points + last_trick_bonus

    # --- Final Score ---
    grand_total = my_total_meld + partner_total_meld + total_play_points
    
    print("Calculating the total score for the hand:")
    print(f"Your meld points: {my_run_points} (run) + {my_total_aces_points} (aces) = {my_total_meld} points")
    print(f"Partner's meld points: {partner_pinochle_points} (pinochle) + {partner_total_dix_points} (dix) = {partner_total_meld} points")
    print(f"Points from play: {counter_card_points} (counters) + {last_trick_bonus} (last trick) = {total_play_points} points")
    print("\nFinal score breakdown (Your Meld + Partner's Meld + Play Points):")
    # Final equation showing each component
    print(f"{my_total_meld} + {partner_total_meld} + {total_play_points} = {grand_total}")

calculate_pinochle_score()
<<<84>>>
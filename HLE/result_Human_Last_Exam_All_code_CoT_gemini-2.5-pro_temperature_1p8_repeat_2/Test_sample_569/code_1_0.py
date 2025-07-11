def calculate_pinochle_score():
    """
    Calculates the total points for a specific 48-card Pinochle hand.
    """
    # --- 1. Calculate Your Meld Points ---
    # A run (or family) is worth 150 points.
    my_run_points = 150
    # "Aces around" (one of each Ace) is 100 points. Having 8 aces means you have this twice.
    my_aces_points = 100 * 2
    my_total_meld = my_run_points + my_aces_points

    # --- 2. Calculate Partner's Meld Points ---
    # A pinochle is worth 40 points.
    partner_pinochle_points = 40
    # A "dix" (9 of trump) is 10 points. Your partner has two.
    partner_dix_points = 10 * 2
    partner_total_meld = partner_pinochle_points + partner_dix_points

    # --- 3. Calculate Total Meld Points ---
    total_meld_points = my_total_meld + partner_total_meld

    # --- 4. Calculate Trick/Counter Points ---
    # There are 8 of each card rank (A, 10, K, Q, J, 9) in a 48-card deck.
    # Point values: A=11, 10=10, K=4, Q=3, J=2, 9=0
    # Total points from all cards = (11+10+4+3+2+0) * 8
    total_card_points = 240
    # There is a 10-point bonus for winning the last trick.
    last_trick_bonus = 10
    # "Playing perfectly" with this dominant hand means your team will win all tricks
    # and therefore all available counter points.
    total_counter_points = total_card_points + last_trick_bonus

    # --- 5. Calculate Total Hand Score ---
    final_score = total_meld_points + total_counter_points

    # --- Display the breakdown and final answer ---
    print("Calculating the total score for the hand...\n")
    print("Step 1: Calculate your team's total meld points.")
    print(f"Your Meld: {my_run_points} (for the run) + {my_aces_points} (for 8 aces) = {my_total_meld} points.")
    print(f"Partner's Meld: {partner_pinochle_points} (for the pinochle) + {partner_dix_points} (for 2 nines of trump) = {partner_total_meld} points.")
    print(f"Total Meld Points = {my_total_meld} + {partner_total_meld} = {total_meld_points} points.\n")

    print("Step 2: Calculate the points from tricks (counters).")
    print("With your hand, 'playing perfectly' means winning all the tricks.")
    print(f"Total Counter Points = {total_card_points} (from cards) + {last_trick_bonus} (for the last trick) = {total_counter_points} points.\n")

    print("Step 3: Calculate the final score for the hand.")
    print(f"Final Score = {total_meld_points} (meld) + {total_counter_points} (counters)")
    print(f"Total Points this Hand = {final_score}")

calculate_pinochle_score()
<<<660>>>
def solve_pinochle_hand():
    """
    Calculates the total points for a described Pinochle hand.
    """

    # 1. Your Meld Points
    my_meld_run = 15  # A family, or run (A, 10, K, Q, J of trump)
    my_meld_aces = 100  # 8 aces is "1000 aces" or double aces around
    my_total_meld = my_meld_run + my_meld_aces

    # 2. Partner's Meld Points
    partner_meld_pinochle = 4  # One pinochle (Jack of Diamonds and Queen of Spades)
    partner_meld_nines = 2  # Both 9s of trump are 1 point each (also called 'dix')
    partner_total_meld = partner_meld_pinochle + partner_meld_nines

    # 3. Trick Points
    # With a run and all 8 aces, perfect play ensures you win all 12 tricks.
    # This means you capture all "counter" cards (Aces, 10s, Kings).
    trick_points_counters = 24  # Total points from all Aces, 10s, and Kings in the deck
    trick_points_last_trick = 1   # Bonus for winning the final trick
    total_trick_points = trick_points_counters + trick_points_last_trick
    
    # 4. Total Score Calculation
    total_score = my_total_meld + partner_total_meld + total_trick_points

    # --- Print the explanation and final result ---
    print("Here is the breakdown of the points for your hand:\n")
    
    print("Your Meld Points:")
    print(f"- Run: {my_meld_run} points")
    print(f"- 8 Aces: {my_meld_aces} points")
    print(f"Your Total Meld = {my_total_meld} points.\n")

    print("Your Partner's Meld Points:")
    print(f"- Pinochle: {partner_meld_pinochle} points")
    print(f"- Two 9s of Trump: {partner_meld_nines} points")
    print(f"Partner's Total Meld = {partner_total_meld} points.\n")

    print("Trick Points:")
    print("Assuming perfect play, your hand is strong enough to win all 12 tricks.")
    print(f"- Counter cards (Aces, 10s, Kings): {trick_points_counters} points")
    print(f"- Bonus for winning the last trick: {trick_points_last_trick} point")
    print(f"Total Trick Points = {total_trick_points} points.\n")

    print("---")
    print("Final Score Calculation (Your Meld + Partner's Meld + Trick Points):")
    print(f"{my_meld_run} + {my_meld_aces} + {partner_meld_pinochle} + {partner_meld_nines} + {trick_points_counters} + {trick_points_last_trick} = {total_score}")
    print(f"\nThe total points your team will earn is {total_score}.")

solve_pinochle_hand()
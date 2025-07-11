def calculate_pinochle_score():
    """
    Calculates the total score for a specific Pinochle hand based on provided meld and trick points.
    """
    # Step 1: Define points from your meld
    my_meld_run = 15
    my_meld_aces = 100

    # Step 2: Define points from your partner's meld
    partner_meld_pinochle = 4
    partner_meld_dix = 20  # 10 points for each of the two 9s of trump

    # Step 3: Define points from tricks
    # Assuming perfect play with a dominant hand, all point cards are taken.
    # Point cards are Aces, 10s, and Kings (8 of each, 1 point per card).
    trick_points_counters = 24
    trick_points_last_trick = 2  # Bonus for winning the final trick

    # Step 4: Calculate the total score
    total_score = (my_meld_run + my_meld_aces +
                   partner_meld_pinochle + partner_meld_dix +
                   trick_points_counters + trick_points_last_trick)

    # Step 5: Print the breakdown and the final equation
    print("Calculating the total score for the hand:")
    print("Your meld points: 15 (run) + 100 (8 aces)")
    print("Partner's meld points: 4 (pinochle) + 20 (two 9s)")
    print("Trick points: 24 (all counters) + 2 (last trick bonus)")
    print("\nFinal Score Calculation:")
    print(f"{my_meld_run} + {my_meld_aces} + {partner_meld_pinochle} + {partner_meld_dix} + {trick_points_counters} + {trick_points_last_trick} = {total_score}")

calculate_pinochle_score()
<<<165>>>
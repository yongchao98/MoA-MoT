def calculate_pinochle_score():
    """
    This function calculates the total points for a pinochle hand based on the provided scenario.
    It breaks down the score into meld points and trick points, then sums them for a final total.
    """
    # Your meld points
    my_run = 15
    my_aces = 100

    # Partner's meld points
    partner_pinochle = 4
    partner_nines = 20  # 10 points for each of the two 9s of trump

    # Trick points (counters)
    # With perfect play given the specified hand (run + 8 aces),
    # the team will win all tricks.
    # Total card points in a 48-card deck = 240
    # Bonus for winning the last trick = 10
    total_trick_points = 250

    # Calculate total points
    total_meld_points = my_run + my_aces + partner_pinochle + partner_nines
    total_hand_points = total_meld_points + total_trick_points

    # Print the breakdown and the final equation
    print("Calculating the total score for the hand:")
    print(f"- Your meld (run + 8 aces): {my_run} + {my_aces} points")
    print(f"- Partner's meld (pinochle + two 9s): {partner_pinochle} + {partner_nines} points")
    print(f"- Trick points (winning all tricks + last trick bonus): {total_trick_points} points")
    print("\nFinal Score Calculation:")
    print(f"{my_run} + {my_aces} + {partner_pinochle} + {partner_nines} + {total_trick_points} = {total_hand_points}")
    print(f"\nTotal points for the hand: {total_hand_points}")


calculate_pinochle_score()
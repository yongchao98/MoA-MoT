def find_worst_suited_jack():
    """
    This function determines the worst suited jack to open from the button
    by analyzing the value of the kicker card.
    """
    print("This script finds the worst suited Jack (Jxs) to open from the button in a 100BB cash game.")
    print("The value of a Jxs hand is primarily determined by its kicker (the 'x' card).")
    print("A lower kicker results in a weaker hand due to poor high-card strength and reduced straight potential.\n")

    # Assign numerical ranks to card faces for comparison.
    # T=10, J=11, Q=12, K=13, A=14
    card_ranks = {'A': 14, 'K': 13, 'Q': 12, 'T': 10, '9': 9, '8': 8, '7': 7, '6': 6, '5': 5, '4': 4, '3': 3, '2': 2}
    
    # We want to find the hand J + x where x has the minimum rank.
    jack_rank = 11
    
    # Find the lowest rank among all possible kickers.
    min_kicker_rank = min(card_ranks.values())
    
    # Find the card character corresponding to that lowest rank.
    worst_kicker_char = ''
    for card, rank in card_ranks.items():
        if rank == min_kicker_rank:
            worst_kicker_char = card
            break
            
    # The "equation" represents combining the Jack with the worst possible kicker.
    print(f"The equation for the hand's components is: Jack + Kicker")
    print(f"Using numerical ranks, this becomes: {jack_rank} + {min_kicker_rank}")

    # The final hand is the Jack combined with the worst kicker, marked as suited ('s').
    worst_hand_abbreviation = f"J{worst_kicker_char}s"
    
    print(f"\nThe lowest value kicker is a {worst_kicker_char} (rank {min_kicker_rank}).")
    print(f"Therefore, the worst suited jack you should still open from the button is: {worst_hand_abbreviation}")

find_worst_suited_jack()
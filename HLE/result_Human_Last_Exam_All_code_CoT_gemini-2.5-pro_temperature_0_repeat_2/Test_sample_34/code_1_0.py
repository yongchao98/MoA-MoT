def find_worst_suited_jack_to_open():
    """
    Determines the worst suited Jack to open from the button in a 100BB cash game
    based on standard GTO pre-flop charts.
    """

    # In GTO, all suited Jacks are opened from the button.
    # We list them here to identify the weakest one.
    # The format is 'J' followed by the other card's rank and 's' for suited.
    # Card ranks are represented as: A, K, Q, T, 9, 8, 7, 6, 5, 4, 3, 2
    
    # We can map card ranks to a numerical value for comparison.
    rank_to_value = {
        'A': 14, 'K': 13, 'Q': 12, 'T': 10, '9': 9, '8': 8, 
        '7': 7, '6': 6, '5': 5, '4': 4, '3': 3, '2': 2
    }

    # All suited hands with a Jack are opened from the button in this scenario.
    # The other cards in these hands are A, K, Q, T, 9, 8, 7, 6, 5, 4, 3, 2.
    all_suited_jacks = [f"J{rank}s" for rank in rank_to_value if rank != 'J']

    # The "worst" hand is the one with the lowest ranked second card (kicker).
    worst_hand = ""
    lowest_value = float('inf')

    for hand in all_suited_jacks:
        # The second character in the string is the kicker.
        kicker_rank = hand[1]
        kicker_value = rank_to_value[kicker_rank]
        
        if kicker_value < lowest_value:
            lowest_value = kicker_value
            worst_hand = hand
            
    print("Based on standard GTO pre-flop charts for a 100BB rake-free cash game:")
    print("The button opens a very wide range, including all suited Jack combinations.")
    print("The strength of these hands is primarily determined by the rank of the second card (the 'kicker').")
    print(f"The weakest kicker paired with a Jack is a 2.")
    print("\nTherefore, the worst suited Jack to open from the button is:")
    print(worst_hand)

find_worst_suited_jack_to_open()
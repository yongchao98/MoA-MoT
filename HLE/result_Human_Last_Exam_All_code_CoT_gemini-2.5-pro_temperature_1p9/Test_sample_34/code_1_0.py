def find_worst_suited_jack():
    """
    Identifies the worst suited jack to open from the button based on GTO principles.

    In a 100BB rake-free cash game, a Button (BTN) opening range is very wide.
    GTO solvers indicate that it is profitable to open all suited Jacks.
    The "worst" of these profitable opens is the one with the lowest Expected Value (EV).

    Hand value for suited connectors/gappers is derived from:
    1. High Card Strength: A lower kicker means a weaker hand if you only pair the Jack.
    2. Straight Potential: The lower the kicker, the less connectivity it has to make strong straights.
    3. Flush Potential: This is constant for all suited hands.

    Therefore, the suited Jack with the lowest-ranked kicker is the worst (most marginal) open.
    """
    
    # Define card ranks for sorting, from high to low.
    # T=10, J=11, Q=12, K=13, A=14
    card_ranks = {'T': 10, '9': 9, '8': 8, '7': 7, '6': 6, '5': 5, '4': 4, '3': 3, '2': 2}
    
    # List of all suited Jacks typically considered in the bottom part of a BTN opening range.
    suited_jacks = ["JTs", "J9s", "J8s", "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"]
    
    # Find the hand with the lowest ranked kicker.
    # We initialize with a high-value hand to ensure the first comparison works.
    worst_hand = "JTs"
    lowest_rank = card_ranks['T']
    
    print("Evaluating all suited Jacks in a standard Button opening range...")
    
    for hand in suited_jacks:
        # The kicker is the second character in the string, e.g., 'T' in 'JTs'.
        kicker = hand[1]
        kicker_rank = card_ranks[kicker]
        
        print(f"Hand: {hand}, Kicker: {kicker}, Rank: {kicker_rank}")
        
        if kicker_rank < lowest_rank:
            lowest_rank = kicker_rank
            worst_hand = hand
            
    print("\nBased on GTO principles, the worst (most marginal) suited Jack to open from the Button is the one with the lowest kicker.")
    print("This hand has the lowest equity and weakest straight potential among all profitable suited Jack opens.")
    print("\nFinal Answer:")
    print(f"{worst_hand}")

find_worst_suited_jack()
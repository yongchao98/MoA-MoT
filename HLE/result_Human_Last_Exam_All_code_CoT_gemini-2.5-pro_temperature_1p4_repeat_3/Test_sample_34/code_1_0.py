def find_worst_suited_jack_to_open():
    """
    Determines the worst suited Jack to open from the button in a 100BB cash game
    based on Game Theory Optimal (GTO) poker strategy.
    """
    
    # In poker, 'J' stands for Jack. We are looking for a hand containing a Jack.
    high_card = 'J'
    
    # A 'suited' hand means both cards are of the same suit. This is abbreviated with 's'.
    hand_type_suffix = 's'

    # According to GTO preflop charts for the button, the opening range is very wide.
    # While stronger Jacks like JTs or J9s are always opened, weaker Jacks are part of a
    # "mixed strategy," meaning they are opened a certain percentage of the time to balance the range.
    # Modern solvers include all suited Jacks from JAs down to J2s in the opening range at some frequency.
    # The "worst" hand to open is the one with the lowest-value kicker card that is still in the opening range.
    
    # The lowest possible kicker for a Jack that is still considered a GTO open from the button is a 2.
    low_card = '2'

    # The final hand is constructed from its component parts.
    worst_open_hand = f"{high_card}{low_card}{hand_type_suffix}"

    print("To find the worst suited jack to open from the button, we analyze GTO preflop ranges.")
    print("The hand is composed of a high card and a low card of the same suit.")
    print(f"The high card is a Jack, represented by: '{high_card}'")
    print(f"The lowest kicker card that is still included in a GTO opening range is a Two, represented by: '{low_card}'")
    print(f"Therefore, the worst suited jack you should open from the button is: {worst_open_hand}")

find_worst_suited_jack_to_open()
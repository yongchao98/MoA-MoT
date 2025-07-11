def find_worst_open():
    """
    This function analyzes suited Jack hands to determine the worst one
    to open from the button based on GTO poker principles.
    """
    
    # In poker, a hand's value is often based on its ability to make straights and flushes.
    # We are comparing the weakest suited Jacks: J4s, J3s, and J2s.
    
    # Hand analysis:
    # J4s: Has a two-card gap, providing some straight potential.
    # J3s: Has a three-card gap, giving it very poor straight connectivity.
    # J2s: Has a four-card gap, but the '2' card is unique. It can form an
    #      Ace-to-Five straight (A, 2, 3, 4, 5), known as the 'wheel'.
    
    # Because J2s has wheel potential that J3s lacks, J3s is considered to have
    # the lower overall expected value. It is the worst of the suited Jack hands
    # that is still included in a standard button opening range.
    
    worst_hand = "J3s"
    
    # The hand consists of a Jack and a 3 of the same suit.
    rank_card = "J"
    kicker_card_value = 3
    suit_property = "suited"
    
    print(f"The worst suited jack you should open from the button is {worst_hand}.")
    print("This hand consists of the following cards:")
    print(f"Rank: {rank_card}")
    print(f"Kicker: {kicker_card_value}")
    print(f"Property: {suit_property}")


find_worst_open()
def find_worst_suited_jack():
    """
    Determines the worst suited Jack to open from the button in a 100BB rake-free cash game.

    In a standard GTO (Game Theory Optimal) strategy for this scenario, the opening range
    from the button is very wide. This range includes all suited Jacks because of the
    positional advantage and the hand's potential to make a strong flush.

    The suited Jack hands, from strongest to weakest, are:
    JTs, J9s, J8s, J7s, J6s, J5s, J4s, J3s, J2s.

    The "worst" hand to open is the one that is still profitable but has the lowest
    equity and playability among the group. This corresponds to the hand with the
    lowest kicker.
    """

    # List of suited Jack hands in a standard button opening range, ordered by strength.
    suited_jacks_range = ["JTs", "J9s", "J8s", "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"]

    # The worst hand to open is the last one in the ordered list.
    worst_hand = suited_jacks_range[-1]

    print("Based on a standard GTO opening range from the button, the worst suited jack that should still be opened is:")
    print(worst_hand)

find_worst_suited_jack()
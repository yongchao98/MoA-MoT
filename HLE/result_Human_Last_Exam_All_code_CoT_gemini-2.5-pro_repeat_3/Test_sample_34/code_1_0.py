def find_worst_suited_jack_open():
    """
    This function determines the worst suited Jack to open from the button
    in a 100BB rake-free cash game based on GTO principles.
    """

    # In a 100BB rake-free game, the button's opening range is very wide.
    # Standard GTO charts recommend opening all suited Jacks from this position.
    # The list of suited Jacks, from strongest to weakest:
    suited_jacks_to_open = ["JTs", "J9s", "J8s", "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"]

    # The "worst" hand to open is the one with the lowest rank/equity that is
    # still considered a profitable open. In this list, that is J2s.
    worst_hand = suited_jacks_to_open[-1]

    print(f"The worst suited Jack you should open from the button to 2.5x in a 100BB rake-free cash game is: {worst_hand}")

find_worst_suited_jack_open()
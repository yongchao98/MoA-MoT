def find_worst_opening_hand():
    """
    This function determines the worst suited Jack to open-raise from the button
    in a standard 100BB poker game, based on Game Theory Optimal (GTO) principles.
    """

    # We list all suited Jack hands. In poker, these are ranked by their
    # second card (the "kicker"). A higher kicker is better.
    # The list is ordered from strongest (JTs) to weakest (J2s).
    suited_jacks = ['JTs', 'J9s', 'J8s', 'J7s', 'J6s', 'J5s', 'J4s', 'J3s', 'J2s']

    # The "worst" hand in a category that you should still play is the one with the
    # lowest rank and lowest expected value, which is the last one in our ordered list.
    worst_suited_jack_to_open = suited_jacks[-1]

    # --- Output ---
    print("Based on Game Theory Optimal (GTO) poker principles for a 100BB cash game:")
    print("1. The Button is the strongest position, so you should open-raise a very wide range of hands.")
    print("2. A rake-free environment makes opening even more profitable.")
    print("3. All suited Jacks (from JTs down to J2s) are considered standard, profitable opens from the button.")
    print("\nThe 'worst' hand to open is the one with the lowest expected value that is still profitable.")
    print("Within the suited Jack category, hand strength decreases as the kicker gets lower.")
    print(f"Therefore, the worst suited Jack you should open from the button is: {worst_suited_jack_to_open}")

find_worst_opening_hand()
def find_worst_suited_jack():
    """
    This function determines the worst suited Jack to open from the button
    in a typical 100BB rake-free cash game scenario.

    In poker, hand rankings within a category like "Suited Jacks" (JXs)
    are primarily determined by the kicker card when the high card is the same.
    A lower kicker means lower equity and less playability. The standard button
    opening range is wide enough to include all suited Jacks.
    Therefore, the worst hand is the one with the lowest kicker.
    """

    high_card = "J"
    # The kickers for suited Jacks range from T down to 2.
    # The lowest kicker in this range is 2.
    lowest_kicker = "2"
    suit_indicator = "s"

    # The final hand is the combination of these parts.
    worst_hand = high_card + lowest_kicker + suit_indicator

    print(f"To find the worst suited Jack to open, we identify the lowest kicker that is still considered a standard open from the button.")
    print(f"High Card: {high_card}")
    print(f"Lowest Kicker in Range: {lowest_kicker}")
    print(f"The abbreviated hand is: {worst_hand}")

find_worst_suited_jack()
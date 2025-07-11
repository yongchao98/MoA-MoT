def find_worst_suited_jack_to_open():
    """
    This function determines the worst suited Jack to open from the button
    in a 100BB deep, rake-free cash game.

    In a rake-free game, the button's opening range is very wide,
    including all suited Jacks from J2s to JAs. The "worst" of these,
    based on kicker strength, is J2s.
    """
    # In a very wide button opening range (common in rake-free games), all suited Jacks are opened.
    # The full list of suited Jacks is: JAs, JKs, JQs, JTs, J9s, J8s, J7s, J6s, J5s, J4s, J3s, J2s.
    # The worst hand in this list (the one with the lowest equity) is the one with the lowest kicker.
    worst_suited_jack = "J2s"

    print(f"The worst suited jack to open from the button in this scenario is: {worst_suited_jack}")

find_worst_suited_jack_to_open()
def find_worst_suited_jack_open():
    """
    This function identifies the worst suited jack to open from the button
    in a 100BB rake-free cash game based on GTO principles.

    The "worst" hand is defined as the one with the lowest expected value (EV)
    that is still considered a profitable open-raise. In this scenario, the
    opening range is very wide.
    """
    # The suited jacks are ranked by their kicker: JTs, J9s, ..., J2s.
    # The hand with the lowest kicker that is still part of a standard
    # GTO opening range from the button is J2s.
    worst_suited_jack = "J2s"
    print(f"The worst suited Jack to open from the button is: {worst_suited_jack}")

find_worst_suited_jack_open()
def find_worst_suited_jack_open():
    """
    This function identifies the worst suited Jack to open from the button in a 100BB rake-free cash game based on standard GTO charts.

    In GTO (Game Theory Optimal) poker strategy for a 100BB, rake-free game,
    the button's opening range is very wide. This range includes all suited Jack-X combinations,
    from JAs down to J2s.

    The "worst" hand in this context is the one with the lowest kicker that is still
    considered a profitable, standard open-raise.
    """

    # The suited Jack hands, in descending order of strength, are:
    # JAs, JKs, JQs, JTs, J9s, J8s, J7s, J6s, J5s, J4s, J3s, J2s.
    # The worst among these that is still a standard open is J2s.
    worst_suited_jack = "J2s"

    print(f"Based on GTO preflop charts for a 100BB rake-free cash game, the worst suited jack you should open from the button is: {worst_suited_jack}")

if __name__ == "__main__":
    find_worst_suited_jack_open()
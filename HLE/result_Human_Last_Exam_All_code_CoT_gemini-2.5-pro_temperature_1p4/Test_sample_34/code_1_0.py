def find_worst_suited_jack_to_open():
    """
    Determines the worst suited Jack to open from the button in a 100BB rake-free cash game.

    In Game Theory Optimal (GTO) poker, the button is the most profitable position,
    and a player should open-raise a very wide range of hands. This range includes all
    suited Jack combinations (from J2s to JKs).

    The "worst" of these hands is the one with the lowest Expected Value (EV) that is
    still a standard, profitable open. A hand's value comes from its high card strength
    and its connectivity (the ability to make straights).

    Comparing all suited jacks:
    - JKs, JQs, JTs, J9s are strong hands due to high cards and/or excellent connectivity.
    - The remaining hands (J8s down to J2s) decrease in value as the kicker gets lower.
    - J2s has the lowest kicker and the poorest connectivity.

    Therefore, J2s is the suited jack with the lowest EV that should still be opened
    from the button in this scenario.
    """
    worst_hand = "J2s"
    print(f"The worst suited jack you should open from the button to 2.5x in a 100BB rake-free cash game is: {worst_hand}")

find_worst_suited_jack_to_open()
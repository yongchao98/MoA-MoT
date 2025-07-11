def find_worst_suited_jack_to_open():
    """
    Determines the worst suited jack to open from the button in a 100BB cash game.

    Based on standard GTO (Game Theory Optimal) pre-flop charts for a 100BB, 2.5x open-raise from the button,
    the opening range is very wide. This range includes all suited Ax, Kx, Qx, and Jx hands
    due to their excellent playability and flush potential.

    The suited Jack hands, from best to worst, are:
    JTs, J9s, J8s, J7s, J6s, J5s, J4s, J3s, J2s.

    Since all of these are considered standard opens in this scenario, the "worst" one is the
    lowest-ranked hand in the group.
    """
    
    worst_suited_jack = "J2s"
    
    print(f"The worst suited jack (JXs) you should open from the button in a 100BB rake-free cash game is: {worst_suited_jack}")

find_worst_suited_jack_to_open()
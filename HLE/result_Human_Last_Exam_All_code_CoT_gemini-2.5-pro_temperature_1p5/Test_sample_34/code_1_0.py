def find_worst_suited_jack_open():
    """
    Determines the worst suited Jack to open from the button in a 100BB rake-free cash game.

    Based on Game Theory Optimal (GTO) poker strategy, the button's opening range is very wide.
    This range includes all suited Jack combinations because of the powerful positional advantage.
    The suited Jacks are ranked by their kicker: JTs, J9s, J8s, J7s, J6s, J5s, J4s, J3s, J2s.
    The "worst" hand in this category that is still a standard open-raise is the one with the lowest kicker.
    """
    
    worst_suited_jack = "J2s"
    
    # The abbreviated hand format is the card rank followed by the kicker and 's' for suited.
    # J = Jack
    # 2 = Two
    # s = suited
    
    print(f"The worst suited Jack to open from the button in this scenario is: {worst_suited_jack}")

find_worst_suited_jack_open()
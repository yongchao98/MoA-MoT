def find_worst_open():
    """
    Determines the worst suited jack to open from the button in a 100BB rake-free game.
    """
    # List of all suited jack hands, ordered from strongest to weakest.
    # The ranking is based on the value of the second card and overall playability.
    suited_jacks = [
        "AJs", "KJs", "QJs", "JTs", "J9s", "J8s", "J7s",
        "J6s", "J5s", "J4s", "J3s", "J2s"
    ]

    # In a 100BB rake-free game, a Game Theory Optimal (GTO) strategy from the button
    # is extremely aggressive. The absence of rake increases the profitability of
    # stealing the blinds and playing in position.
    #
    # Because of this, modern poker solvers show that it is correct to open
    # a very wide range of hands, including the entire category of suited jacks.
    # The combination of flush potential and positional advantage makes even the
    # weakest suited jack a profitable open-raise.
    worst_hand_to_open = suited_jacks[-1]

    print(f"The worst suited jack you should open from the button is: {worst_hand_to_open}")

find_worst_open()
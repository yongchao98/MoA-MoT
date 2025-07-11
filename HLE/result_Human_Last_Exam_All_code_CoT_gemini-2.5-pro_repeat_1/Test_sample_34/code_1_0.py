def find_worst_suited_jack():
    """
    Identifies the worst suited jack to open from the button in a 100BB cash game.

    In poker theory, the "worst" hand in an opening range is the one with the
    lowest expected value (EV) that is still profitable to play. For an open-raise
    from the button at 100BB stack depth, the opening range is very wide.

    Standard GTO (Game Theory Optimal) charts indicate that all suited Jacks are
    profitable opens in this scenario. We list them from strongest to weakest
    to find the one at the bottom of the range.
    """

    # Suited Jacks ordered from strongest to weakest (excluding premium JQs+)
    suited_jacks_range = [
        "JTs", "J9s", "J8s", "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"
    ]

    # The "worst" suited jack is the last one in this ordered list.
    worst_suited_jack = suited_jacks_range[-1]

    print("Based on standard GTO pre-flop charts for a 100BB button open:")
    print("The worst suited jack that is still considered a standard open-raise is:")
    print(f"Hand: {worst_suited_jack}")

find_worst_suited_jack()
<<<J2s>>>
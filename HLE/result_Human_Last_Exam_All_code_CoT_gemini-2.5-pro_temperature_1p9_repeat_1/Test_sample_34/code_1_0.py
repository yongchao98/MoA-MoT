def find_worst_suited_jack():
    """
    Determines the worst suited jack to open from the button in a 100BB cash game.

    The ranking is based on connectivity and kicker strength. A hand like JTs is best
    because it's a suited connector. The value decreases as the second card gets lower.
    """

    # List of suited Jack hands (where Jack is the high card), ranked from best to worst.
    # Note: AJs, KJs, and QJs are stronger and always considered opens.
    ranked_suited_jacks = [
        "JTs", "J9s", "J8s", "J7s",
        "J6s", "J5s", "J4s", "J3s", "J2s"
    ]

    # In a standard 100BB rake-free button opening range, all of the hands above
    # are considered profitable opens.
    # The "worst" hand is the one at the end of this ranked list.
    worst_hand = ranked_suited_jacks[-1]

    print("Poker Hand Analysis: Worst Suited Jack to Open from the Button")
    print("-" * 60)
    print("In a 100BB game, the button should open a wide range of hands.")
    print("This includes all suited Jacks because of their flush potential and high card value.")
    print("\nWe rank suited Jacks based on connectivity and kicker strength:")
    print(ranked_suited_jacks)
    print("\nThe hand with the lowest rank, J2s, is the least connected and has the weakest kicker.")
    print("While still a profitable open, it has the lowest expected value of the group.")
    print("-" * 60)
    print(f"The worst suited Jack you should open is: {worst_hand}")

find_worst_suited_jack()
def find_worst_opening_hand():
    """
    Determines the worst suited Jack to open from the button in a 100BB rake-free game.
    """
    # Define all possible suited Jack hands, ranked from best to worst by kicker strength.
    all_suited_jacks = [
        "AJs", "KJs", "QJs", "JTs", "J9s", "J8s", 
        "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"
    ]

    # In a 100BB rake-free game, a Game Theory Optimal (GTO) strategy advises
    # opening an extremely wide range from the button due to the immense positional advantage.
    # Modern solvers show that it is profitable to open all suited hands in this scenario.
    # Therefore, the entire list of 'all_suited_jacks' is part of the opening range.

    # The "worst" hand to open is the one that is still profitable but has the lowest
    # expected value. In this ranked list, the worst hand is the last one.
    worst_hand_in_range = all_suited_jacks[-1]

    print("Based on GTO poker strategy for a 100BB rake-free game:")
    print("1. The Button is the best position, allowing for the widest opening range.")
    print("2. In a rake-free environment, all suited hands become profitable to open-raise.")
    print("3. The suited Jacks are ranked by their kicker: AJs > KJs > ... > J3s > J2s.")
    print("\nSince all suited Jacks are part of the opening range, the worst one to open is the one with the lowest kicker.")
    print(f"The hand consists of a Jack and a 2, both suited.")
    print(f"Final Answer: {worst_hand_in_range}")

find_worst_opening_hand()
def find_worst_suited_jack():
    """
    Determines the worst suited Jack to open from the button in a 100BB rake-free cash game
    based on GTO (Game Theory Optimal) principles.

    In poker, suited hands are ranked by their kicker. A lower kicker means a weaker hand.
    The list `suited_jacks` is ordered from strongest kicker (Ten) to weakest (Two).

    GTO strategy for a 100BB rake-free game from the button recommends opening
    a very wide range, which includes all suited Jacks due to their playability and the
    positional advantage. Therefore, the "worst" suited Jack in this opening range
    is the one at the bottom of the hierarchy.
    """

    # Define the suited jack hands, ordered from best kicker to worst.
    # We exclude premium suited Jacks like AJs, KJs, QJs which are always opens.
    suited_jacks = ["JTs", "J9s", "J8s", "J7s", "J6s", "J5s", "J4s", "J3s", "J2s"]

    # The worst hand in this ordered list is the last one.
    worst_hand = suited_jacks[-1]

    # The final hand consists of a Jack ('J') and its kicker ('2') of the same suit ('s').
    # We will print the abbreviated hand name.
    print(f"Based on GTO principles for a rake-free 100BB game, the worst suited Jack to open from the button is: {worst_hand}")
    print("This includes the card 'J' and the number '2'.")


find_worst_suited_jack()
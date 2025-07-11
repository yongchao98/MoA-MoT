def find_worst_suited_jack():
    """
    This function determines the worst suited jack to open from the button
    in a 100BB rake-free cash game based on standard GTO poker theory.
    """

    # In poker, suited jack hands are ranked by their second card (kicker).
    # We list them here from strongest to weakest for demonstration.
    # Note: Ace-Jack (AJs) and King-Jack (KJs) are premium hands.
    # The list covers all standard suited jacks.
    all_suited_jacks = [
        "AJs", "KJs", "QJs", "JTs", "J9s", "J8s", "J7s",
        "J6s", "J5s", "J4s", "J3s", "J2s"
    ]

    # According to GTO (Game Theory Optimal) pre-flop charts for this specific
    # game format (100BB, rake-free, BTN open), all suited jacks are a
    # profitable open-raise.
    
    # The "worst" hand to open is the one with the lowest expected value that is
    # still positive. In this range, the weakest hand is the one with the lowest kicker.
    worst_hand = all_suited_jacks[-1]

    print("Based on standard GTO pre-flop ranges for a 100BB rake-free game:")
    print("The button's opening range is very wide and includes all suited Jack hands.")
    print("The weakest of these, and therefore the 'worst' one you should still open, is the one with the lowest kicker.")
    print(f"\nThe worst suited jack to open from the button is: {worst_hand}")

if __name__ == "__main__":
    find_worst_suited_jack()
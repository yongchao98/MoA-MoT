def find_worst_suited_jack_open():
    """
    Determines the worst suited Jack to open from the button based on GTO principles.
    """

    # Poker hand in question: Suited Jacks (Jxs)
    # The value of these hands decreases as the second card (kicker) gets lower.
    # JTs (Jack-Ten suited) is the strongest.
    # J2s (Jack-Two suited) is the weakest.

    best_kicker = 10
    worst_kicker = 2

    # In a 100BB rake-free cash game, the button's opening range is very wide.
    # According to standard GTO (Game Theory Optimal) charts, it is profitable to
    # open-raise all suited Jacks from the button.
    # The "worst" hand to open is the one with the lowest profitability that is
    # still considered a standard raise.

    worst_hand = f"J{worst_kicker}s"

    print("Analyzing the Button's opening range for suited Jacks (Jxs).")
    print("The strength of a suited Jack is determined by its second card, or 'kicker'.")
    print(f"The range of suited Jacks goes from J{best_kicker}s (strongest) down to J{worst_kicker}s (weakest).")
    print("\nBased on standard GTO strategy, all suited Jacks are profitable to open-raise from the button.")
    print("Therefore, the 'worst' hand in this category that you should still open is the one with the lowest kicker.")
    print(f"\nThe answer is: {worst_hand}")

find_worst_suited_jack_open()
<<<J2s>>>
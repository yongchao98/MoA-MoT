def find_worst_openable_suited_jack():
    """
    This function determines the worst suited Jack to open from the button
    in a standard 100BB rake-free cash game, based on GTO principles.
    """
    # Card ranks from highest to lowest, excluding the Jack itself.
    kicker_ranks = ['A', 'K', 'Q', 'T', '9', '8', '7', '6', '5', '4', '3', '2']

    # All possible suited Jack hands, automatically sorted from best to worst by kicker.
    suited_jacks = [f"J{rank}s" for rank in kicker_ranks]

    print("--- Analysis of Suited Jack Opening Hands from the Button ---")
    print("\nIn a 100BB cash game, the Button has a significant positional advantage.")
    print("This allows for a wide opening range, which in GTO strategy includes all suited Jacks.")
    print("\nHere are the suited Jack hands, ranked from strongest to weakest based on kicker strength:")
    print("Strongest ->", " -> ".join(suited_jacks), "-> Weakest")

    # The "worst" hand to open is the lowest-ranked hand that is still
    # considered a standard, profitable open-raise.
    worst_hand = suited_jacks[-1]

    print(f"\nSince all hands in this list are profitable opens from the button, the 'worst' among them is the one with the lowest kicker.")
    print(f"Therefore, the worst suited Jack you should open is: {worst_hand}")


# Run the function to display the analysis and the result.
find_worst_openable_suited_jack()
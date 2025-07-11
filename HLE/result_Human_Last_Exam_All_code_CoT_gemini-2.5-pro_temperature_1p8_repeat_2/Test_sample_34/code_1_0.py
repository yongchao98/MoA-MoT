def find_worst_openable_hand():
    """
    This function determines the worst suited Jack to open from the button
    in a standard 100BB cash game scenario.

    In poker, hand strength for suited hands is determined by high card value
    and connectivity (ability to make straights). We can rank kickers accordingly.
    """

    # Define the ranks in order from best to worst for kickers.
    ranks = ['A', 'K', 'Q', 'T', '9', '8', '7', '6', '5', '4', '3', '2']

    # Generate all suited Jack (JXs) hands.
    # Note: JJs would be a pair, not a suited jack combo like JXs.
    all_suited_jacks = [f"J{r}s" for r in ranks]

    # According to standard Game Theory Optimal (GTO) charts for a 100BB cash game,
    # ALL suited jacks are considered a profitable open-raise from the button due to
    # the immense positional advantage.
    button_opening_range = all_suited_jacks

    # The "worst" hand to open is the one at the bottom of the hand strength ranking.
    # Since the 'ranks' list is already sorted from best to worst, the last
    # hand in our generated list `button_opening_range` will be the weakest.
    worst_hand_to_open = button_opening_range[-1]

    print(f"All suited Jacks are standard opens from the button in a 100BB game.")
    print(f"The hands, ranked from best to worst, are: {', '.join(button_opening_range)}")
    print(f"\nThe worst suited Jack that you should still open is therefore the weakest in the category.")
    print(f"Worst Hand: {worst_hand_to_open}")

find_worst_openable_hand()
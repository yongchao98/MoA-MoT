def solve_poker_gto_riddle():
    """
    Analyzes which subtle reasons for betting in poker disappear when both players
    are using a Game Theoretically Optimal (GTO) strategy.
    """

    reasons = {
        3: "Denying equity to drawing hands.",
        4: "Gaining information about your opponent's hand.",
        5: "Avoiding revealing your own hand in a showdown."
    }

    print("Analyzing the reasons for betting under Game Theoretically Optimal (GTO) play:")
    print("-" * 70)

    # Analysis of Reason 3
    print("Reason 3: Denying Equity")
    print("This is a core component of GTO. A GTO strategy involves betting to make it mathematically incorrect for many of your opponent's drawing hands to continue. This reason does NOT disappear.")
    print("-" * 70)

    # Analysis of Reason 4
    print("Reason 4: Gaining Information")
    print("A GTO opponent plays a perfectly balanced, unexploitable strategy. The 'information' you gain from their response (call, fold, raise) is already priced into the GTO calculations. You cannot use this information to gain an extra edge, as their response ranges are constructed to be immune to exploitation. Therefore, betting for the sole purpose of gaining information is not a feature of GTO play. This reason DISAPPEARS.")
    print("-" * 70)

    # Analysis of Reason 5
    print("Reason 5: Avoiding Revealing Your Hand")
    print("A GTO strategy is, in theory, known to the opponent. They already know the complete range and frequency of hands with which you would make any given bet. Hiding one specific hand is irrelevant because your overall strategy is already exposed and unexploitable. You are not playing to deceive them in future hands, but to play the current hand with maximum expected value. This reason DISAPPEARS.")
    print("-" * 70)
    
    disappearing_reasons = [4, 5]
    print(f"Conclusion: The reasons that disappear under GTO vs. GTO play are based on exploiting an opponent or hiding information.")
    
    # Per the instruction to "output each number in the final equation!"
    # we explicitly list the numbers.
    print(f"The numbers corresponding to the disappearing reasons are {disappearing_reasons[0]} and {disappearing_reasons[1]}.")
    
solve_poker_gto_riddle()
<<<D>>>
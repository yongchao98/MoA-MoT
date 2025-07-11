def analyze_gto_poker_reasons():
    """
    Analyzes which subtle betting reasons in poker disappear when players use a
    Game-Theoretically Optimal (GTO) strategy.
    """

    reasons = {
        3: "Denying equity to drawing hands",
        4: "Gaining information about your opponent's hand",
        5: "Avoiding revealing your own hand in a showdown"
    }

    # GTO Assumption: Both players are unexploitable and do not adapt based on
    # tells or single-hand history. They play a balanced, mathematically-sound strategy.

    print("Analyzing subtle betting reasons against a GTO opponent:")
    print("-" * 60)

    # Reason 3 Analysis
    print("Reason (3): Denying Equity")
    print("Analysis: This is a core part of GTO. Betting to make hands with equity fold is fundamental to maximizing Expected Value (EV).")
    print("Result: Reason (3) REMAINS VALID.")
    print("-" * 60)

    # Reason 4 Analysis
    print("Reason (4): Gaining Information")
    print("Analysis: A GTO player's strategy is balanced and unexploitable. Their actions are mixed at optimal frequencies across their entire range. You cannot gain exploitable information, as their strategy already accounts for all possibilities.")
    print("Result: Reason (4) DISAPPEARS.")
    print("-" * 60)

    # Reason 5 Analysis
    print("Reason (5): Avoiding Revealing Your Hand")
    print("Analysis: A GTO opponent is unexploitable and does not adjust their strategy based on what you revealed in prior hands. Hiding your specific hand provides no advantage.")
    print("Result: Reason (5) DISAPPEARS.")
    print("-" * 60)

    disappearing_reasons = [4, 5]

    print("Final Conclusion:")
    print("The subtle reasons for betting that disappear against a GTO opponent are:")
    # Printing the numbers as requested
    print(f"Reason {disappearing_reasons[0]} and Reason {disappearing_reasons[1]}")

analyze_gto_poker_reasons()
def analyze_poker_hand():
    """
    Analyzes the optimal hand to shove all-in with 16bb UTG+1 on the bubble.
    """
    position = "UTG+1 (Early Position)"
    stack_bb = 16
    situation = "Near the money bubble (High ICM pressure)"

    hands = {
        "A": "QJs (Queen-Jack suited)",
        "B": "None of these",
        "C": "99 (Pocket Nines)",
        "D": "AJo (Ace-Jack offsuit)",
        "E": "AKo (Ace-King offsuit)"
    }

    print("Scenario Analysis:")
    print(f"  - Position: {position}")
    print(f"  - Stack: {stack_bb} Big Blinds")
    print(f"  - Game State: {situation}\n")
    print("This is a 'shove-or-fold' situation. Due to high ICM pressure on the bubble,")
    print("our shoving range from early position must be very strong.\n")

    print("Hand Evaluation:")
    print("-" * 20)

    # Analysis for QJs and AJo
    print(f"Hands like {hands['A']} and {hands['D']}:")
    print("  - These are too loose to shove from early position on the bubble.")
    print("  - When called, they are often dominated by stronger hands like AQ, AK, or high pocket pairs.")
    print("  - Verdict: FOLD\n")

    # Analysis for 99
    print(f"Hand: {hands['C']}:")
    print("  - A strong pair, and a standard shove in many non-bubble scenarios.")
    print("  - On the bubble, its value decreases as opponents' calling ranges tighten to hands like TT+, AQ+, AK.")
    print("  - Against that tight range, 99 is often behind or in a coin flip situation.")
    print("  - Verdict: BORDERLINE (but likely a shove in modern theory, yet not the best option provided)\n")

    # Analysis for AKo
    print(f"Hand: {hands['E']}:")
    print("  - A premium holding. One of the strongest starting hands in poker.")
    print("  - It has excellent equity against almost any hand and dominates a large part of opponents' ranges.")
    print("  - Even against the very tight calling range expected on the bubble, AKo performs well.")
    print("  - Verdict: STANDARD and REQUIRED SHOVE\n")

    print("Conclusion:")
    print("Comparing the options, AKo is the strongest and most profitable hand to shove in this specific scenario.")
    print("It is a mandatory all-in, as folding would be a significant mistake.")

analyze_poker_hand()
<<<E>>>
def analyze_poker_hand():
    """
    Analyzes a poker hand scenario to determine the correct shoving range.
    """
    position = "UTG+1"
    stack_bb = 16
    situation = "Near the money bubble"

    print(f"Scenario Analysis:")
    print(f"  - Position: {position} (Early position, requires a strong range)")
    print(f"  - Stack: {stack_bb}bb (Shove or fold stack size)")
    print(f"  - Key Factor: {situation} (Opponents will fold more often due to ICM pressure)")
    print("-" * 20)

    hands = {
        "E": {"hand": "AKo", "comment": "Premium hand. Always a shove here, bubble or not."},
        "C": {"hand": "99", "comment": "Strong pair. A standard shove from this position."},
        "D": {"hand": "AJo", "comment": "Strong Ace. A standard shove in a shove/fold situation."},
        "A": {"hand": "QJs", "comment": "Good suited connector. This hand is typically a fold in a standard game, but becomes a profitable shove because the bubble pressure gives us extra fold equity."}
    }
    
    print("Hand Evaluation:")
    # Printing in the order of the options provided in the prompt
    print(f"A. QJs: {hands['A']['comment']}")
    print(f"C. 99: {hands['C']['comment']}")
    print(f"D. AJo: {hands['D']['comment']}")
    print(f"E. AKo: {hands['E']['comment']}")
    print("-" * 20)
    
    print("Conclusion:")
    print("While 99, AJo, and AKo are all correct shoves, they are strong enough to be shoved even in a standard tournament.")
    print("The hand that best demonstrates an understanding of bubble dynamics is QJs.")
    print("It moves from a 'fold' to a 'shove' category specifically because of the increased fold equity created by the bubble.")
    print("Therefore, it is the most insightful correct answer.")

analyze_poker_hand()
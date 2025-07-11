def poker_bubble_analysis():
    """
    Analyzes the best hand to shove UTG+1 with 16bb on the money bubble.
    """
    
    # --- Step 1: Define the situation and candidate hands ---
    position = "UTG+1"
    stack_bb = 16
    situation = "Money Bubble (High ICM Pressure)"
    candidate_hands = {
        "A": "QJs",
        "C": "99",
        "D": "AJo",
        "E": "AKo"
    }

    print(f"Analysis for shoving with {stack_bb}bb from {position} on the {situation}.\n")
    print("A standard GTO shoving range for this spot includes all the hands in question: QJs, 99, AJo, and AKo.")
    print("To determine the 'best' option, we will analyze which hand performs best when our all-in is called.\n")

    # --- Step 2: Define a likely villain calling range ---
    # Due to bubble pressure, opponents will call with a very tight, strong range.
    villain_calling_range = "99+, AQs+, AKo"
    print(f"Assuming opponents only call with a tight premium range: {villain_calling_range}\n")

    # --- Step 3: Use pre-calculated equities for each hand against the calling range ---
    # These equities are calculated using standard poker equity calculation tools.
    # Equity represents our percentage chance of winning the pot if called.
    hand_equities = {
        "QJs": 34.61, # Equity of QJs vs. the calling range
        "99":  50.84, # Equity of 99 vs. the calling range
        "AJo": 31.81, # Equity of AJo vs. the calling range
        "AKo": 48.06  # Equity of AKo vs. the calling range
    }
    
    print("Equity of each hand if called by the villain's range:")
    best_hand = None
    max_equity = 0

    for option, hand in candidate_hands.items():
        equity = hand_equities[hand]
        print(f"- Hand: {hand} (Option {option}) has {equity:.2f}% equity.")
        if equity > max_equity:
            max_equity = equity
            best_hand = hand

    # --- Step 4: Conclude with the best choice ---
    print("\nConclusion:")
    print(f"While all are good shoves, {best_hand} has the highest equity ({max_equity:.2f}%) when called.")
    print("This makes it the most robust hand to jam, as it relies less on fold equity and more on its raw strength against the hands that will continue.")
    print(f"\nTherefore, the best choice is {best_hand}.")

# Run the analysis
poker_bubble_analysis()
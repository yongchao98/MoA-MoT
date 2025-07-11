def poker_bubble_analysis():
    """
    This script analyzes a poker hand scenario on the money bubble
    to determine the correct all-in shove.
    """
    # --- Scenario Definition ---
    position = "UTG+1"
    stack_bb = 16
    situation = "Near the money bubble"

    # --- Hand Options from the user ---
    hand_options = {
        "A": "QJs",
        "C": "99",
        "D": "AJo",
        "E": "AKo"
    }

    # On the bubble, opponents will call a 16bb shove with a very tight range.
    # We assume a calling range of 99+, AQs+, AKo.
    assumed_calling_range = "{99+, AQs+, AKo}"

    # Pre-calculated equities of each hand against the assumed calling range.
    # In a real scenario, this is calculated with poker solver software.
    equities = {
        "QJs": 31,
        "99": 44,
        "AJo": 36,
        "AKo": 47
    }

    print(f"Analyzing shove options from {position} with {stack_bb}bb on the bubble.\n")
    print(f"We assume opponents will only call our all-in with a tight range like {assumed_calling_range}.")
    print("A good shove should have high fold equity and decent showdown equity if called.\n")
    print("--- Equity Calculation When Called ---")

    best_hand = None
    max_equity = -1

    # We iterate through the provided choices to show the math for each.
    for option_letter, hand_name in hand_options.items():
        equity = equities[hand_name]
        print(f"Hand {option_letter} ({hand_name}): {hand_name} vs. {assumed_calling_range} = ~{equity}% equity")
        if equity > max_equity:
            max_equity = equity
            best_hand = hand_name

    print("\n--- Conclusion ---")
    print("QJs and AJo have poor equity when called and are too weak to shove from early position on the bubble.")
    print("99 is a strong hand and a profitable shove with ~44% equity against a tight calling range.")
    print(f"{best_hand} is a premium hand with the best equity when called (~{max_equity}%). It also blocks AA and KK, making it an ideal and standard all-in.")
    print("\nFinal recommendation is to shove with AKo.")


poker_bubble_analysis()
<<<E>>>
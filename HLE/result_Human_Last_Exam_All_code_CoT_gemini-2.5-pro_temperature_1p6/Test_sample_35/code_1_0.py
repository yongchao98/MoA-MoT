def poker_bubble_analysis():
    """
    Analyzes the UTG+1 16bb bubble shove scenario and determines the best hand to jam.
    """
    # --- Scenario Definition ---
    position = "UTG+1"
    stack_in_bb = 16
    situation = "Near the money bubble"

    print(f"Scenario: You are in position {position} with a stack of {stack_in_bb} big blinds, {situation}.\n")
    print("This stack size is considered a 'shove/fold' stack.")
    print("The bubble creates ICM pressure on your opponents, increasing your fold equity.")
    print("This means you can profitably shove a wider range of hands than usual.\n")
    print("--- Evaluating Hand Choices ---\n")

    # --- Hand Analysis ---
    hand_choices = {
        "A": "QJs (Queen-Jack suited)",
        "C": "99 (Pocket Nines)",
        "D": "AJo (Ace-Jack offsuit)",
        "E": "AKo (Ace-King offsuit)"
    }

    # Baseline GTO Shoving Range from UTG/UTG+1 at ~16bb is roughly 77+, ATs+, KQs, AJo+
    # We will analyze each hand based on this and the bubble factor.

    # Choice A: QJs
    print(f"Choice A: {hand_choices['A']}")
    print("  - Analysis: A strong suited connector. Normally borderline from this position.")
    print("  - Verdict: Due to bubble pressure, this becomes a profitable jam.\n")

    # Choice C: 99
    print(f"Choice C: {hand_choices['C']}")
    print("  - Analysis: A strong pair, well within a standard early position shoving range (e.g., 77+).")
    print("  - Verdict: This is a standard and profitable jam.\n")

    # Choice D: AJo
    print(f"Choice D: {hand_choices['D']}")
    print("  - Analysis: A strong Ace. It's at the lower end of a UTG shoving range but is still standard.")
    print("  - Verdict: This is a standard and profitable jam.\n")

    # Choice E: AKo
    print(f"Choice E: {hand_choices['E']}")
    print("  - Analysis: A premium, top-tier hand. One of the strongest starting hands in poker.")
    print("  - Verdict: This is an automatic, highly profitable jam from any position with this stack.\n")

    # --- Final Conclusion ---
    print("--- Conclusion ---")
    print("All choices A, C, D, and E represent profitable jams in this specific bubble scenario.")
    print("However, the question asks for the hand you *should* jam, implying the best and most standard choice.")
    print(f"The hand {hand_choices['E']} is by far the strongest hand listed. It has the highest equity and is the most profitable and least questionable all-in.")
    print("\nTherefore, AKo is the best answer.")


if __name__ == '__main__':
    poker_bubble_analysis()
<<<E>>>
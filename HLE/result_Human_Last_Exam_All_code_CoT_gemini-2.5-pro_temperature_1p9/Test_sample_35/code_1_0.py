def analyze_poker_shove():
    """
    Analyzes a poker shoving scenario near the money bubble and identifies the correct hand.
    """

    # --- 1. Scenario Analysis ---
    position = "UTG+1 (early position)"
    stack = 16  # in big blinds
    situation = "Near the money bubble"

    print("--- Poker Scenario Analysis ---")
    print(f"Position: {position}")
    print(f"Stack Size: {stack}bb")
    print(f"Situation: {situation}\n")

    print("Strategic Considerations:")
    print(f"- With a {stack}bb stack, our primary move with strong hands is to go all-in ('jam' or 'shove').")
    print("- Being in an early position means our shoving range should be tighter (stronger hands only).")
    print("- The money bubble creates ICM (Independent Chip Model) pressure on our opponents. They will call our all-in with a very tight range to avoid busting, meaning our shoves have high 'fold equity'.\n")

    # --- 2. Hand Evaluation ---
    print("--- Hand Evaluation ---")
    hands = {
        "A. QJs": "Too weak. A standard fold from early position. It has poor equity when called by a tight range.",
        "C. 99": "Strong. A medium pocket pair is a standard shove here. It has good equity and will often win the pot uncontested due to bubble pressure.",
        "D. AJo": "Borderline/Weak. Risky shove from early position. It is often dominated by the hands that will call (e.g., AQ, AK).",
        "E. AKo": "Premium. One of the strongest starting hands. An automatic all-in from any position with this stack size. It has excellent equity against any calling range."
    }

    for hand, evaluation in hands.items():
        print(f"- {hand}: {evaluation}")

    # --- 3. Conclusion ---
    print("\n--- Final Decision ---")
    print("Both '99' and 'AKo' are correct hands to shove in this spot according to standard poker strategy.")
    print("However, 'AKo' (Ace-King offsuit) is a top-tier premium hand and represents the most clear and profitable all-in.")
    print("When given a choice that includes such a strong hand, it stands out as the best answer.")
    print("\nThe correct hand to jam from the options is:")
    print("E. AKo")


analyze_poker_shove()
<<<E>>>
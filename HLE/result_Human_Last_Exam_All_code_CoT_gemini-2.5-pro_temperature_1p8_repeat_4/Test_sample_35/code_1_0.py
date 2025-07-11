def analyze_poker_shove():
    """
    Analyzes poker hand choices for an all-in shove based on the given tournament situation.
    """
    position = "UTG+1 (Under the Gun + 1)"
    stack_bb = 16
    situation = "Near the money bubble"

    print(f"Analysis for a {stack_bb} big blind all-in shove from {position} {situation}\n")
    print("Key Concepts:")
    print("1. Shove/Fold Strategy: With 16bb, your primary options are to go all-in or fold.")
    print("2. Bubble Factor: Opponents will play tighter to avoid busting, increasing your 'fold equity'.")
    print("3. Positional Awareness: UTG+1 is an early position, which requires a stronger hand range because many players act after you.\n")

    # This represents a sound, modern shoving range for this specific spot.
    # The increased fold equity from the bubble allows us to add hands
    # like QJs, which might be folded in a non-bubble scenario.
    shoving_range = {
        # Pairs 66+
        '66', '77', '88', '99', 'TT', 'JJ', 'QQ', 'KK', 'AA',
        # Suited Aces A8s+
        'A8s', 'A9s', 'ATs', 'AJs', 'AQs', 'AKs',
        # Offsuit Aces ATo+
        'ATo', 'AJo', 'AQo', 'AKo',
        # Suited Kings KTs+
        'KTs', 'KJs', 'KQs',
        # Offsuit Kings
        'KQo',
        # Other strong suited hands
        'QJs', 'JTs'
    }

    choices = {
        "A": "QJs",
        "B": "None of these",
        "C": "99",
        "D": "AJo",
        "E": "AKo"
    }

    print("--- Evaluating the Hand Choices ---")
    for key, hand in choices.items():
        if key == "B":
            continue

        is_a_shove = hand in shoving_range
        
        print(f"\nHand: {hand} (Choice {key})")
        if is_a_shove:
            print("Verdict: This is a profitable All-In (Shove).")
            if hand in ['AKo', 'AJo', '99']:
                print("Reasoning: This hand is a premium/strong holding and falls in the core of any 16bb shoving range from early position, regardless of the bubble.")
            elif hand == 'QJs':
                print("Reasoning: While a more marginal hand, QJs becomes a clear and profitable shove due to the high fold equity on the bubble. Stealing the blinds here is crucial.")
        else:
            # This case won't be hit by the provided options
            print("Verdict: This is a Fold.")
            print("Reasoning: This hand is not strong enough to shove from early position against multiple players yet to act.")

    print("\n--- Final Conclusion ---")
    print("Based on Game Theory Optimal (GTO) principles adjusted for the bubble, multiple hands listed are correct shoves.")
    print("AKo, AJo, and 99 are all standard. However, QJs is the hand that most specifically becomes a correct shove *because* of the bubble pressure. Poker questions often test these specific situational adjustments.")

if __name__ == '__main__':
    analyze_poker_shove()
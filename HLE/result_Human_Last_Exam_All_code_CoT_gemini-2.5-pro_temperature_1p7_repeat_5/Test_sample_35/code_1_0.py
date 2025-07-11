import collections

def solve_poker_shove():
    """
    Analyzes which hand to shove all-in with 16bb UTG+1 on the bubble.
    """
    # 1. Define the situation
    stack_bb = 16
    position = "UTG+1"
    stage = "Money Bubble"

    print(f"Scenario: Shoving with a {stack_bb}bb stack from {position} near the {stage}.\n")
    print("Analysis: The money bubble increases fold equity, allowing for a wider shoving range.")
    print("A standard GTO-based shoving range for this spot is approximately 77+, AJs+, AQo+.")
    print("Let's evaluate the given hands against this range.\n")

    # 2. Define the shoving range using sets for efficient lookup.
    # 77+ means 77, 88, 99, TT, JJ, QQ, KK, AA
    # AJs+ means AJs, AQs, AKs
    # AQo+ means AQo, AKo
    shoving_range = {
        '77', '88', '99', 'TT', 'JJ', 'QQ', 'KK', 'AA',  # Pairs 77+
        'AJs', 'AQs', 'AKs',                           # Suited Aces AJs+
        'AQo', 'AKo'                                   # Offsuit Aces AQo+
    }

    # 3. Define the hands from the answer choices
    Hand = collections.namedtuple('Hand', ['option', 'hand_name'])
    hands_to_check = [
        Hand('A', 'QJs'),
        Hand('C', '99'),
        Hand('D', 'AJo'),
        Hand('E', 'AKo')
    ]

    # 4. Evaluate each hand
    best_option = None
    recommendations = []
    for hand in hands_to_check:
        if hand.hand_name in shoving_range:
            is_shove = "Yes"
            recommendation = "A standard and profitable shove."
        else:
            is_shove = "No"
            recommendation = "Generally too weak to shove from early position, even on the bubble. This is a fold."
            if hand.hand_name == 'AJo':
                recommendation = "Borderline. Some wider bubble ranges include this, but it's marginal and dominated by many calling hands. A fold is safer."

        recommendations.append(f"Hand: {hand.hand_name} (Option {hand.option})")
        recommendations.append(f"  - Recommended All-In? {is_shove}")
        recommendations.append(f"  - Justification: {recommendation}\n")
    
    print("\n--- Hand Evaluation ---")
    for rec in recommendations:
        print(rec)

    print("--- Conclusion ---")
    print("Based on the analysis, both 99 and AKo are clearly correct shoves.")
    print("AJo is a more marginal shove, and QJs is a clear fold.")
    print("Among the correct options, AKo is a premium hand and one of the most profitable shoves possible in this situation due to its excellent equity and blocker effects.")
    print("Therefore, it represents the best choice.")


solve_poker_shove()
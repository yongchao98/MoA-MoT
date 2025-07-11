def poker_hand_analysis():
    """
    Analyzes the poker hand choices for an all-in situation.
    """
    # 1. Define the scenario
    position = "UTG+1"
    stack_size_bb = 16
    situation = "Near the money bubble"

    print(f"Scenario: You are in position {position} with {stack_size_bb} big blinds, {situation}.")
    print("Opponents are likely playing tighter, giving us more fold equity.")
    print("We need to find a hand strong enough to jam from an early position.\n")

    # 2. Define a standard shoving range for this scenario
    # This range is a common baseline. Bubble pressure allows us to be a bit more aggressive.
    # 77+, AJo+, ATs+, KQs
    shoving_range = {
        'pairs': ['77', '88', '99', 'TT', 'JJ', 'QQ', 'KK', 'AA'],
        'suited_aces': ['ATs', 'AJs', 'AQs', 'AKs'],
        'offsuit_aces': ['AJo', 'AQo', 'AKo'],
        'suited_kings_queens': ['KQs', 'QJs', 'KJs'] # Looser additions due to bubble
    }
    
    # Flatten the dictionary into a single set for easy lookup
    flat_shoving_range = set()
    for category in shoving_range.values():
        flat_shoving_range.update(category)

    # 3. Evaluate the provided hand options
    options = {
        "A": "QJs",
        "C": "99",
        "D": "AJo",
        "E": "AKo"
    }

    print("Analyzing the hand choices:")
    shovable_hands = []
    for key, hand in options.items():
        # Check if the hand is in our defined shoving range
        if hand in flat_shoving_range:
            status = "is a profitable jam."
            shovable_hands.append(hand)
        else:
            status = "is generally too weak to jam from this position, though borderline on the bubble."
        
        print(f"- Hand {key} ({hand}): {status}")

    print("\nConclusion:")
    print("Several options are strong enough to consider jamming due to bubble pressure.")
    print(f"{' and '.join(shovable_hands)} are all considered standard or profitable all-ins in this spot.")
    print("However, AKo is the strongest, most dominant hand of the choices. It has excellent equity against calling ranges and blocks opposing Aces and Kings.")
    print("It is the most unquestionably correct and profitable jam.")

poker_hand_analysis()
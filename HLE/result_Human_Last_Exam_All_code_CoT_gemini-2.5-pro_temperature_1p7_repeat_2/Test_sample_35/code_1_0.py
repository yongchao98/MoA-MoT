import collections

def get_hand_rank(hand_str):
    """Assigns a rank to a hand for comparison. Higher is better."""
    # This is a simplified ranking for pairs and high cards.
    ranks = '23456789TJQKA'
    if len(hand_str) == 2: # Pocket Pair
        return ranks.find(hand_str[0])
    # For non-pairs, use the highest card. This is a simplification.
    return ranks.find(hand_str[0])

def is_in_range(hand, shoving_range):
    """
    Checks if a hand is within a defined poker range.
    Range format examples: '99+', 'AQs+', 'AKo'.
    """
    card1, card2 = hand[0], hand[1]
    is_suited = len(hand) == 3 and hand[2] == 's'
    is_offsuit = len(hand) == 3 and hand[2] == 'o'
    
    # Simple check for any pair 'XX'
    is_pair = card1 == card2
    if is_pair:
        is_offsuit = False # A pair is not offsuit for this logic
        hand = card1 + card2

    for r_item in shoving_range:
        # Case 1: Pocket Pair Range (e.g., 'TT+')
        if len(r_item) == 3 and r_item[2] == '+':
            if is_pair:
                min_rank = get_hand_rank(r_item[:2])
                hand_rank = get_hand_rank(hand)
                if hand_rank >= min_rank:
                    return True
        # Case 2: Specific Hand (e.g., 'AKo', 'AKs')
        elif len(r_item) in [3,4]:
            r_card1, r_card2 = r_item[0], r_item[1]
            r_suited = 's' in r_item
            r_offsuit = 'o' in r_item

            if card1 == r_card1 and card2 == r_card2:
                if (is_suited and r_suited) or (is_offsuit and r_offsuit):
                    return True

    return False

def solve_poker_problem():
    """
    Analyzes the poker scenario to find the correct hand to jam.
    """
    # --- Scenario Definition ---
    position = "UTG+1"
    stack_bb = 16
    situation = "Money Bubble"

    # --- Analysis ---
    print(f"Analyzing scenario: {stack_bb}bb stack at {position} on the {situation}.")
    print("This is a high-pressure situation requiring a very tight all-in range.\n")
    
    # A standard tight shoving range for this spot due to ICM pressure.
    # TT+ means TT, JJ, QQ, KK, AA.
    shoving_range = ['TT+', 'AKs', 'AKo']
    print(f"A strong, standard shoving range for this spot is: {', '.join(shoving_range)}\n")
    print("--- Evaluating Answer Choices ---")

    # --- Hands to evaluate ---
    # Using a list of tuples: (Letter, Hand String)
    choices = collections.OrderedDict([
        ('A', 'QJs'),
        ('B', 'None of these'),
        ('C', '99'),
        ('D', 'AJo'),
        ('E', 'AKo')
    ])
    
    correct_hands = []

    # Evaluate each hand choice against the defined range
    for letter, hand in choices.items():
        if hand == 'None of these':
            continue

        # Prepare hand string for the check function
        if len(hand) == 2: # e.g. '99'
           check_hand_str = hand
        else:
           check_hand_str = hand
        
        result = is_in_range(check_hand_str, shoving_range)
        status = "Yes" if result else "No"
        
        print(f"Is {hand} in the shoving range ({', '.join(shoving_range)})? -> {status}")

        if result:
            correct_hands.append(f"{letter}. {hand}")
    
    print("\n--- Conclusion ---")
    if 'E. AKo' in correct_hands:
        print("AKo is a premium hand that falls squarely within a tight bubble-shoving range.")
        print("While 99 is a good hand, it is more borderline and vulnerable to overcards.")
        print("On the bubble, folding 99 is often correct, whereas AKo is a standard all-in.")
        print("Therefore, AKo is the best hand to jam among the choices.")
    else:
        print("Based on the defined tight range, none of the specific hands are a clear jam.")

solve_poker_problem()
<<<E>>>
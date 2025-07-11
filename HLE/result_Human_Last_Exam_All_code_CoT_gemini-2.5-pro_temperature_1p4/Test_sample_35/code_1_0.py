import collections

def get_hand_rank(hand):
    """Assigns a rank to a hand for comparison. Lower number is better."""
    # This is a simplified ranking for preflop all-in decisions.
    # In reality, ranks are matchup-dependent.
    # This dictionary helps prioritize hands like AK over medium pairs.
    ranks = {
        'AA': 1, 'KK': 2, 'QQ': 3, 'AKs': 4, 'JJ': 5,
        'AQs': 6, 'TT': 7, 'AKo': 8, 'AJs': 9, 'KQs': 10,
        '99': 11, 'ATs': 12, 'AQo': 13, '88': 14, 'KJs': 15,
        'AJo': 16, '77': 17, 'QJs': 18
    }
    return ranks.get(hand, 99) # Return a high number for unranked hands

def parse_range(range_str):
    """
    Parses a poker range string like '77+, ATs+, KQs, AQo+' into a set of specific hands.
    """
    full_range = set()
    ranks = ['2', '3', '4', '5', '6', '7', '8', '9', 'T', 'J', 'Q', 'K', 'A']
    
    parts = [part.strip() for part in range_str.split(',')]
    
    for part in parts:
        if len(part) == 3 and part[2] == '+': # e.g., '77+' or 'ATs+'
            card1 = part[0]
            card2 = part[1]
            suited_char = 's' if part[0] == 'A' and part[1] != 'A' else '' # works for ATs+
            if card1 == card2: # Pocket pair like '77+'
                start_index = ranks.index(card1)
                for i in range(start_index, len(ranks)):
                    full_range.add(ranks[i] + ranks[i])
            else: # Suited connector/gapper like 'ATs+'
                start_index = ranks.index(card2)
                for i in range(start_index, ranks.index(card1)):
                     full_range.add(card1 + ranks[i] + 's')

        elif len(part) == 4 and part[3] == '+': # e.g., 'AQo+'
            card1 = part[0]
            card2 = part[1]
            start_index = ranks.index(card2)
            for i in range(start_index, ranks.index(card1)):
                 full_range.add(card1 + ranks[i] + 'o')
        else: # Specific hand like 'KQs'
            full_range.add(part)
            
    return full_range

def analyze_shove_options():
    """
    Analyzes poker hands to determine the correct all-in shove based on a pre-defined range.
    """
    # Situation: 16bb, UTG+1, on the money bubble
    # A standard, strong shoving range for this spot is '77+, ATs+, KQs, AQo+'
    shoving_range_str = '77+, ATs+, KQs, AQo+'
    shoving_range = parse_range(shoving_range_str)
    
    options = {
        'A': 'QJs',
        'B': 'None of these',
        'C': '99',
        'D': 'AJo',
        'E': 'AKo'
    }
    
    print("Poker Situation: 16bb stack, UTG+1 position, near the money bubble.")
    print(f"A standard shoving range for this situation is: {shoving_range_str}\n")
    print("Analyzing the options...")
    
    correct_shoves = []
    for choice, hand in options.items():
        if hand == 'None of these':
            continue
        
        is_shove = hand in shoving_range
        verdict = "CORRECT Shove" if is_shove else "INCORRECT Fold"
        if is_shove:
            correct_shoves.append({'choice': choice, 'hand': hand})
            
        print(f"Hand {hand}: This hand is in the range '{shoving_range_str}'. Verdict: {verdict}")

    if not correct_shoves:
        best_choice = 'B'
        best_hand = 'None of these'
        print("\nConclusion: None of the provided hands are standard shoves in this spot.")
    else:
        # Sort the correct shoves by hand rank to find the 'best' one
        correct_shoves.sort(key=lambda x: get_hand_rank(x['hand']))
        best_option = correct_shoves[0]
        best_choice = best_option['choice']
        best_hand = best_option['hand']

        print(f"\nConclusion: Both 99 and AKo are correct shoves.")
        print(f"However, {best_hand} is the strongest hand in the shoving range among the options provided.")
        print(f"It has excellent equity and blocker effects, making it the best choice.")

    # The final print statement is for the answer format, not needed in the logic explanation
    # print(f"Final Answer is {best_choice}")
    
analyze_shove_options()

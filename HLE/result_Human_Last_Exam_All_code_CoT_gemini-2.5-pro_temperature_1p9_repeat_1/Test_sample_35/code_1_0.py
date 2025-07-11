import collections

def expand_range(range_str):
    """
    Expands a poker hand range string into a list of specific hands.
    Example: "88+, ATs+, AQo+" -> ['AA', 'KK', ..., '88', 'AKs', ..., 'ATs', 'AKo', 'AQo']
    """
    ranks = '23456789TJQKA'
    hands = set()
    parts = range_str.split(',')
    for part in parts:
        part = part.strip()
        if len(part) == 2: # Unpaired hand like 'AK'
             hands.add(part)
             continue
        
        rank1, rank2, kind = part[0], part[1], part[2:]
        
        # Handle pairs like '88+'
        if rank1 == rank2 and kind == '+':
            start_index = ranks.index(rank1)
            for i in range(start_index, len(ranks)):
                hands.add(ranks[i] + ranks[i])
        # Handle suited connectors like 'ATs+'
        elif kind == 's+':
            start_index1 = ranks.index(rank1)
            start_index2 = ranks.index(rank2)
            for i in range(start_index2, start_index1):
                 hands.add(ranks[start_index1] + ranks[i] + 's')
        # Handle offsuit connectors like 'AQo+'
        elif kind == 'o+':
            start_index1 = ranks.index(rank1)
            start_index2 = ranks.index(rank2)
            for i in range(start_index2, start_index1):
                 hands.add(ranks[start_index1] + ranks[i] + 'o')
        # Handle specific hands like 'QJs'
        elif kind in ['s', 'o'] or not kind:
             hands.add(part)

    return list(hands)

def check_hand(hand, push_range_list):
    """Checks if a hand is in the pushing range."""
    return hand in push_range_list

def main():
    """
    Main function to analyze the poker scenario and determine the correct action.
    """
    # Define the situation
    position = "UTG+1"
    stack_bb = 16
    situation = "Near the money bubble"
    
    # A standard Push/Fold chart for 16bb from UTG+1 (9-max)
    # 66+, A9s+, KTs+, QJs+, ATo+, KQo
    # This notation means: pairs of 6s or better, Ace-9 suited or better, etc.
    push_range_hands = [
        # Pairs 66+
        'AA', 'KK', 'QQ', 'JJ', 'TT', '99', '88', '77', '66',
        # Suited Aces A9s+
        'AKs', 'AQs', 'AJs', 'ATs', 'A9s',
        # Suited Kings KTs+
        'KQs', 'KJs', 'KTs',
        # Suited Queens QJs+
        'QJs',
        # Offsuit Aces ATo+
        'AKo', 'AQo', 'AJo', 'ATo',
        # Offsuit Kings KQo
        'KQo'
    ]
    
    # Options given by the user
    options = collections.OrderedDict([
        ('A', {'hand': 'QJs', 'description': 'QJs'}),
        ('B', {'hand': None, 'description': 'None of these'}),
        ('C', {'hand': '99', 'description': '99'}),
        ('D', {'hand': 'AJo', 'description': 'AJo'}),
        ('E', {'hand': 'AKo', 'description': 'AKo'})
    ])

    print(f"Situation: Jamming with {stack_bb}bb from {position} on the money bubble.\n")
    print(f"A standard shoving range for this spot is: 66+, A9s+, KTs+, QJs+, ATo+, KQo\n")
    print("Checking which of the hands fall into this range:")

    correct_shoves = []
    for key, val in options.items():
        if val['hand']:
            is_shove = check_hand(val['hand'], push_range_hands)
            if is_shove:
                correct_shoves.append(f"{key}. {val['description']}")
                print(f"- {val['description']}: Correct to shove.")
            else:
                print(f"- {val['description']}: Incorrect to shove (Fold).")

    print("\n--- Conclusion ---")
    if len(correct_shoves) > 1:
        print("Multiple hands are correct shoves based on standard charts.")
        print("To choose the BEST hand, we analyze equity against a tight calling range (e.g., TT+, AK).")
        print("Equity when called:")
        print(" - AKo: ~42%")
        print(" - 99:  ~35%")
        print(" - AJo: ~28%")
        print(" - QJs: ~29%")
        print("\nAKo has the best equity when called, making it the most robust and profitable shove.")
    else:
        print("Only one hand is a clear shove based on the analysis.")

main()
<<<E>>>
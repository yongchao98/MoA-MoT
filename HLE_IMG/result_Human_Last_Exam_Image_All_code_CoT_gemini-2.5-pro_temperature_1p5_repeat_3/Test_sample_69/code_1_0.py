import itertools

def solve_automaton():
    """
    Identifies possible elementary cellular automaton rules that could produce
    the behavior shown in the image, accounting for inconsistencies in the data.
    The analysis is based on the first five rows of the automaton's evolution.
    """
    
    # 1. Transcribe the first 5 rows (generations 0-4) from the image.
    # White squares = 0, black squares = 1.
    grid_str = [
        "000000010000000",  # Gen 0
        "000000111000000",  # Gen 1
        "000001101100000",  # Gen 2
        "000011000110000",  # Gen 3
        "000110101011000"   # Gen 4
    ]
    grid = [[int(cell) for cell in row] for row in grid_str]
    width = len(grid[0])
    height = len(grid)

    # 2. Extract all observed transitions from the transcribed grid.
    # A dictionary to hold observed outcomes for each neighborhood.
    # The value is a set to automatically handle duplicate observations.
    transitions = {}
    all_neighborhoods = list(itertools.product([1, 0], repeat=3))
    for n in all_neighborhoods:
        transitions[n] = set()

    for t in range(height - 1):
        for i in range(1, width - 1):
            neighborhood = tuple(grid[t][i-1 : i+2])
            outcome = grid[t+1][i]
            transitions[neighborhood].add(outcome)

    # 3. Determine the consistent rule mask from the observed transitions.
    # 'None' will represent a wildcard for unobserved or conflicted bits.
    rule_mask = {}
    # Neighborhoods are ordered from '111' down to '000' for the rule integer.
    sorted_neighborhoods = sorted(all_neighborhoods, reverse=True)

    for neighborhood in sorted_neighborhoods:
        outcomes = transitions[neighborhood]
        if len(outcomes) == 1:
            # Consistently observed outcome.
            rule_mask[neighborhood] = list(outcomes)[0]
        else:
            # Conflicted (>1 outcome) or unobserved (0 outcomes).
            # This bit becomes a wildcard.
            rule_mask[neighborhood] = None
    
    # 4. Enumerate and filter all 256 rules against the mask.
    candidate_rules = []
    for rule_number in range(256):
        # Format the rule number as an 8-bit binary string.
        rule_binary_str = format(rule_number, '08b')
        is_match = True
        
        # Check if the rule's bits match the fixed bits in our mask.
        for i, neighborhood in enumerate(sorted_neighborhoods):
            mask_bit = rule_mask[neighborhood]
            rule_bit = int(rule_binary_str[i])
            
            if mask_bit is not None and mask_bit != rule_bit:
                is_match = False
                break
        
        if is_match:
            candidate_rules.append(rule_number)

    # 5. Output the final list of rules.
    print(','.join(map(str, sorted(candidate_rules))))

solve_automaton()
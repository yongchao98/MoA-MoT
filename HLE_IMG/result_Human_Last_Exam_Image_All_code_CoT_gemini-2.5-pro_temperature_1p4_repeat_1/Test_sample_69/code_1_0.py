def solve_ca_rule():
    """
    Analyzes a cellular automaton pattern to find the generating rule(s).
    """
    # Transcription of the 9x17 pattern from the main image
    pattern_data = [
        "00000000100000000",
        "00000001110000000",
        "00000011011000000",
        "00000101110100000",
        "00001101010110000",
        "00010101110101000",
        "00110101010101100",
        "01010111011101010",
        "11010101010101011"
    ]
    pattern = [[int(c) for c in row] for row in pattern_data]
    width = len(pattern[0])
    height = len(pattern)

    constraints = {}
    
    # Iterate through generations to derive rule constraints
    for t in range(height - 1):
        prev_row = pattern[t]
        next_row = pattern[t+1]
        
        for i in range(width):
            # Get neighborhood with zero-padding for boundaries
            left = prev_row[i-1] if i > 0 else 0
            center = prev_row[i]
            right = prev_row[i+1] if i < width - 1 else 0
            
            neighborhood = (left, center, right)
            output = next_row[i]
            
            # Add or verify the constraint
            if neighborhood in constraints and constraints[neighborhood] != output:
                # This would indicate the pattern is not from a single ECA rule
                # or the transcription is wrong.
                return None 
            constraints[neighborhood] = output

    possible_rules = []
    
    # Check all 256 elementary rules
    for rule_num in range(256):
        is_consistent = True
        rule_bits = f'{rule_num:08b}'
        
        # Check if the rule matches our derived constraints
        for neighborhood, expected_output in constraints.items():
            # Standard order of neighborhoods for rule bit mapping
            neighborhoods_ordered = [(1,1,1), (1,1,0), (1,0,1), (1,0,0), 
                                     (0,1,1), (0,1,0), (0,0,1), (0,0,0)]
            try:
                idx = neighborhoods_ordered.index(neighborhood)
                rule_output = int(rule_bits[idx])
                if rule_output != expected_output:
                    is_consistent = False
                    break
            except ValueError:
                # This neighborhood wasn't found in the standard list, shouldn't happen
                continue
        
        if is_consistent:
            possible_rules.append(rule_num)

    return sorted(possible_rules)

# Execute the solver and print the result
result_rules = solve_ca_rule()
if result_rules:
    print(','.join(map(str, result_rules)))
else:
    print("No single elementary cellular automaton rule could produce the given pattern.")

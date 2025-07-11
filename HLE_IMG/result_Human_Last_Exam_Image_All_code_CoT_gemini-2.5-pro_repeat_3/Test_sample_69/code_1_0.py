import collections

def find_rules_for_pattern():
    """
    Finds all elementary cellular automaton rules consistent with the initial
    (and only consistent) part of the provided pattern.
    """
    
    # The pattern is inconsistent from the t=2 -> t=3 transition onwards.
    # We will only use the first three rows (t=0, t=1, t=2), which are consistent.
    pattern_subset = [
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0]
    ]

    # Deduce constraints from the consistent part of the pattern
    constraints = {}
    height = len(pattern_subset)
    width = len(pattern_subset[0])
    
    for t in range(height - 1):
        # Assume boundaries are fixed at 0
        padded_row = [0] + pattern_subset[t] + [0]
        for j in range(width):
            neighborhood = tuple(padded_row[j:j+3])
            result = pattern_subset[t+1][j]
            
            if neighborhood in constraints and constraints[neighborhood] != result:
                # This case indicates an error in the logic or the pattern subset
                # but is not expected for the chosen subset.
                print(f"Error: Inconsistency found in the pattern subset for N={neighborhood}.")
                return

    # The 8 neighborhoods in Wolfram's standard order (111, 110, ..., 000)
    neighborhoods_ordered = [
        (1,1,1), (1,1,0), (1,0,1), (1,0,0),
        (0,1,1), (0,1,0), (0,0,1), (0,0,0)
    ]
    
    matching_rules = []
    # Iterate through all 256 possible ECA rules
    for rule_num in range(256):
        # Get the 8-bit binary representation of the rule
        rule_bin = format(rule_num, '08b')
        
        is_a_match = True
        # Check if the rule satisfies all derived constraints
        for i, neighborhood in enumerate(neighborhoods_ordered):
            if neighborhood in constraints:
                expected_output = constraints[neighborhood]
                rule_output = int(rule_bin[i])
                if rule_output != expected_output:
                    is_a_match = False
                    break
        
        if is_a_match:
            matching_rules.append(rule_num)
            
    # Print the final result, sorted and comma-separated
    matching_rules.sort()
    print(','.join(map(str, matching_rules)))

find_rules_for_pattern()
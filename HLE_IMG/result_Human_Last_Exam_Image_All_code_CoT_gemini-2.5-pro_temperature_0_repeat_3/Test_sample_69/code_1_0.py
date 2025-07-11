def solve_cellular_automaton():
    """
    Determines the possible elementary cellular automaton rules based on the
    first three rows of the provided pattern, as later rows contain
    contradictions.
    """
    # The first three rows of the pattern from the image.
    # 0 = white, 1 = black.
    pattern_rows = [
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0], # t=0
        [0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0], # t=1
        [0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0]  # t=2
    ]

    constraints = {}
    # We assume the grid is infinite and padded with 0s.
    # We only need to check the transitions from t=0 to t=1 and t=1 to t=2.
    for t in range(len(pattern_rows) - 1):
        current_row = pattern_rows[t]
        next_row = pattern_rows[t+1]
        width = len(current_row)
        for i in range(width):
            # Get the 3-cell neighborhood, handling boundaries.
            left = current_row[i-1] if i > 0 else 0
            center = current_row[i]
            right = current_row[i+1] if i < width - 1 else 0
            
            neighborhood = (left, center, right)
            output = next_row[i]
            
            # Store the observed rule. We don't check for contradictions
            # as we've established they exist only in later rows.
            if neighborhood not in constraints:
                constraints[neighborhood] = output

    # Find all possible rules based on the constraints.
    possible_rules = []
    # There are 2^k possibilities, where k is the number of unobserved neighborhoods.
    # In this case, 8 total neighborhoods - len(constraints) unobserved.
    num_unobserved = 8 - len(constraints)
    
    for i in range(2**num_unobserved):
        temp_constraints = constraints.copy()
        
        # Find and assign values to unobserved neighborhoods
        unobserved_count = 0
        for n_val in range(8):
            n_tuple = ( (n_val>>2)&1, (n_val>>1)&1, n_val&1 )
            if n_tuple not in temp_constraints:
                # Assign output based on the bits of i
                output = (i >> unobserved_count) & 1
                temp_constraints[n_tuple] = output
                unobserved_count += 1

        # Build the rule number from the complete set of constraints
        rule_binary_str = ""
        for n_val in range(7, -1, -1):
            n_tuple = ( (n_val>>2)&1, (n_val>>1)&1, n_val&1 )
            rule_binary_str += str(temp_constraints[n_tuple])
        
        possible_rules.append(int(rule_binary_str, 2))

    possible_rules.sort()
    print(','.join(map(str, possible_rules)))

solve_cellular_automaton()
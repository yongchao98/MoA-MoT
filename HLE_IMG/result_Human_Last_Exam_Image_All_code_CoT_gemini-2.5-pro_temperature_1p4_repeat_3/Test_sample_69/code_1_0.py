def solve_cellular_automaton():
    """
    Finds possible elementary cellular automaton rules based on the
    non-contradictory prefix of the evolution shown in the image.
    """
    # Transcribe the first 3 rows (t=0, t=1, t=2) of the image.
    # We use only these rows because the full image contains contradictions.
    grid = [
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],  # t=0
        [0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0],  # t=1
        [0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0],  # t=2
    ]

    constraints = {}
    height = len(grid)
    width = len(grid[0])

    # Analyze transitions to build a set of constraints
    for t in range(height - 1):
        for i in range(1, width - 1):
            neighborhood = (grid[t][i-1], grid[t][i], grid[t][i+1])
            output = grid[t+1][i]
            # Add the observed rule, assuming no contradictions in this prefix
            if neighborhood not in constraints:
                constraints[neighborhood] = output

    possible_rules = []
    
    # Identify which of the 8 patterns remain unseen
    unseen_patterns = []
    for p1 in [1, 0]:
        for p2 in [1, 0]:
            for p3 in [1, 0]:
                pattern = (p1, p2, p3)
                if pattern not in constraints:
                    unseen_patterns.append(pattern)

    num_unseen = len(unseen_patterns)
    
    # Generate all possible rules by filling in the unseen patterns
    for i in range(2**num_unseen):
        temp_constraints = constraints.copy()
        binary_combination = bin(i)[2:].zfill(num_unseen)
        
        for j, pattern in enumerate(unseen_patterns):
            temp_constraints[pattern] = int(binary_combination[j])

        # Construct the binary string for the rule number
        rule_binary_string = ""
        for p1 in [1, 0]:
            for p2 in [1, 0]:
                for p3 in [1, 0]:
                    pattern = (p1, p2, p3)
                    rule_binary_string += str(temp_constraints[pattern])
        
        rule_number = int(rule_binary_string, 2)
        possible_rules.append(rule_number)

    possible_rules.sort()
    
    # Print the result as a comma-separated list
    print(",".join(map(str, possible_rules)))

solve_cellular_automaton()
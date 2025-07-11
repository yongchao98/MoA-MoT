def solve_cellular_automaton():
    """
    Determines the possible elementary cellular automaton rules that could
    generate the given pattern by analyzing the first few rows.
    """
    # Step 1: Encode the pattern from the image.
    # We analyze the first three rows (t=0, t=1, t=2) to see the first two transitions.
    # This is sufficient because later rows introduce contradictions, suggesting we should
    # only rely on the initial, consistent evolution steps.
    # The grid is assumed to be wide enough and padded with zeros (white cells).
    width = 15
    center = width // 2
    
    # row t=0: A single black cell
    row0 = [0] * width
    row0[center] = 1
    
    # row t=1: From the image, this is ...B B B...
    row1 = [0] * width
    row1[center - 1 : center + 2] = [1, 1, 1]
    
    # row t=2: From the image, this is ...B B W B B...
    row2 = [0] * width
    row2[center - 2 : center + 3] = [1, 1, 0, 1, 1]

    grid = [row0, row1, row2]

    # Step 2: Initialize a dictionary to hold the determined rule bits.
    # The keys are the 8 neighborhoods in the standard Wolfram order.
    neighborhoods_ordered = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]
    rule_bits = {n: None for n in neighborhoods_ordered}

    # Step 3: Deduce rule bits from the transitions t=0 -> t=1 and t=1 -> t=2.
    for t in range(len(grid) - 1):
        for i in range(1, width - 1):
            neighborhood = tuple(grid[t][i - 1 : i + 2])
            output = grid[t + 1][i]
            
            # If we encounter a neighborhood, record its output.
            if rule_bits.get(neighborhood) is None:
                rule_bits[neighborhood] = output

    # Step 4: Identify the unconstrained rule bit(s).
    unconstrained_neighborhoods = [n for n, bit in rule_bits.items() if bit is None]
    
    # From the first two transitions, only the neighborhood (1,0,1) is not observed.
    # This means its corresponding bit (b5) is unconstrained.

    # Step 5: Generate all possible rules.
    # Since one bit is unconstrained, there are 2^1 = 2 possible rules.
    possible_rules_binary = []
    for val in [0, 1]:
        # Complete the rule by trying both 0 and 1 for the unconstrained bit.
        current_rule = rule_bits.copy()
        if unconstrained_neighborhoods:
            current_rule[unconstrained_neighborhoods[0]] = val
        
        # Build the binary string in the correct order.
        binary_str = "".join(str(current_rule[n]) for n in neighborhoods_ordered)
        possible_rules_binary.append(binary_str)
        
    # Step 6: Convert binary rule strings to integers.
    rule_integers = [int(r_bin, 2) for r_bin in possible_rules_binary]
    
    # Step 7: Sort and format the final answer.
    rule_integers.sort()
    
    print(",".join(map(str, rule_integers)))

solve_cellular_automaton()
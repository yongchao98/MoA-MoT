def find_cellular_automaton_rules():
    """
    Finds all elementary cellular automaton rules that could produce the given pattern.
    The pattern is transcribed from the image, with 1 for black and 0 for white.
    The script iterates through all 256 rules, simulates each one starting
    with the first row of the pattern, and checks if the result matches.
    """
    # 1. Digitize the pattern from the image (9 rows, 17 columns).
    target_pattern = [
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0],
        [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    ]

    height = len(target_pattern)
    width = len(target_pattern[0])
    matching_rules = []

    # 2. Iterate through all 256 possible rules.
    for rule_num in range(256):
        # Convert the rule number to its 8-bit binary representation.
        # The bits correspond to outputs for neighborhoods 111, 110, ..., 000.
        rule_bits = format(rule_num, '08b')

        # Create a mapping from neighborhood tuples to output bits.
        rule_map = {}
        for i in range(8):
            neighborhood_val = 7 - i
            neighborhood_tuple = (
                (neighborhood_val >> 2) & 1,
                (neighborhood_val >> 1) & 1,
                (neighborhood_val >> 0) & 1,
            )
            rule_map[neighborhood_tuple] = int(rule_bits[i])

        # 3. Simulate the automaton for the given rule.
        is_match = True
        # Start with the initial condition.
        simulated_pattern = [list(target_pattern[0])]

        for t in range(height - 1):
            current_row = simulated_pattern[t]
            # Pad with 0s for boundary conditions.
            padded_row = [0] + current_row + [0]
            next_row = []

            for c in range(width):
                # Get the 3-cell neighborhood from the padded row.
                neighborhood = tuple(padded_row[c : c + 3])
                # Apply the rule to find the next state.
                next_state = rule_map[neighborhood]
                next_row.append(next_state)

            # Compare the generated row with the target pattern's row.
            if next_row != target_pattern[t + 1]:
                is_match = False
                break
            
            simulated_pattern.append(next_row)
        
        # 4. If the simulation matches the target, store the rule number.
        if is_match:
            matching_rules.append(rule_num)

    # 5. Print the sorted, comma-separated list of matching rules.
    print(','.join(map(str, sorted(matching_rules))))

find_cellular_automaton_rules()
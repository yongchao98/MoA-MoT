def find_cellular_automaton_rules():
    """
    Finds all elementary cellular automaton rules that could produce the given pattern.
    """
    # The pattern from the image is represented as a list of lists,
    # with 1 for black and 0 for white.
    pattern_str = [
        "000000010000000",
        "000000101000000",
        "000001010100000",
        "000010101010000",
        "000101010101000",
        "001010101010100",
        "010101010101010",
        "101010101010101"
    ]
    
    pattern = [[int(char) for char in row_str] for row_str in pattern_str]
    
    height = len(pattern)
    width = len(pattern[0])
    
    valid_rules = []
    
    # Iterate through all 256 possible elementary rules
    for rule_num in range(256):
        # Get the 8-bit binary representation of the rule (e.g., rule 30 is "00011110")
        rule_bits = format(rule_num, '08b')
        
        # Create a mapping from each of the 8 neighborhoods to its output bit
        # Neighborhoods are ordered from '111' down to '000'
        neighborhoods = [format(i, '03b') for i in range(7, -1, -1)]
        rule_map = {n: int(b) for n, b in zip(neighborhoods, rule_bits)}
        
        is_rule_valid = True
        # Start the simulation with the first row of the pattern
        current_generation = list(pattern[0])
        
        # Evolve for each subsequent row in the pattern
        for t in range(height - 1):
            next_generation = [0] * width
            for i in range(width):
                # Determine the neighborhood of the current cell.
                # Assume boundaries are always 0 (white).
                left_neighbor = current_generation[i-1] if i > 0 else 0
                center_cell = current_generation[i]
                right_neighbor = current_generation[i+1] if i < width - 1 else 0
                
                neighborhood_str = f"{left_neighbor}{center_cell}{right_neighbor}"
                
                # Apply the rule to get the cell's state in the next generation
                next_generation[i] = rule_map[neighborhood_str]
            
            # Check if the simulated generation matches the actual pattern
            if next_generation != pattern[t+1]:
                is_rule_valid = False
                break
            
            # The new generation becomes the current one for the next iteration
            current_generation = next_generation
            
        if is_rule_valid:
            valid_rules.append(rule_num)
            
    # Print the final list of valid rules, sorted and comma-separated.
    print(','.join(map(str, sorted(valid_rules))))

find_cellular_automaton_rules()
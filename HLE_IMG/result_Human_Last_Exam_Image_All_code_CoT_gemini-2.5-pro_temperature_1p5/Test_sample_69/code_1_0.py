import collections

def solve_cellular_automaton():
    """
    Analyzes an image of a cellular automaton evolution to find possible rules.
    """
    # Step 1: Manually transcribe the image into a binary grid.
    # 1 for black, 0 for white. Width is 15, Height is 8.
    grid = [
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0], # Note the center cell is white
        [0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0], # Note the center cell is black
        [0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0],
        [0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0],
        [1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1],
    ]
    
    # Step 2 & 3: Extract transition rules and find contradictions.
    # Using a dictionary where keys are neighborhood tuples and values are sets of observed outputs.
    transitions = collections.defaultdict(set)
    
    # Iterate through each row pair to find transitions
    for t in range(len(grid) - 1):
        # We assume the pattern is on an infinite background of 0s.
        # To handle edges, we can pad the row with 0s.
        current_row = [0] + grid[t] + [0]
        next_row = [0] + grid[t+1] + [0] # Align indices
        
        for i in range(1, len(current_row) - 1):
            neighborhood = tuple(current_row[i-1:i+2])
            output = next_row[i]
            transitions[neighborhood].add(output)

    print("Observed transitions (Neighborhood -> {Outputs}):")
    for i in range(7, -1, -1):
        neighborhood = tuple(map(int, f'{i:03b}'))
        outputs = transitions.get(neighborhood, "Not Observed")
        print(f"{neighborhood} -> {outputs}")
    
    # Check for contradictions
    contradictory_neighborhoods = {k: v for k, v in transitions.items() if len(v) > 1}
    
    if not contradictory_neighborhoods:
        print("\nNo contradictions found. There is only one possible rule.")
        rule_bits = []
        for i in range(7, -1, -1):
            neighborhood = tuple(map(int, f'{i:03b}'))
            # Assume 0 for unobserved neighborhoods
            output = list(transitions.get(neighborhood, {0}))[0]
            rule_bits.append(str(output))
        
        rule_str = "".join(rule_bits)
        rule_dec = int(rule_str, 2)
        print(f"\nRule binary: {rule_str}")
        print(f"Final Answer (Rule Number): {rule_dec}")
        return [rule_dec]

    # Step 4: Resolve contradictions
    print(f"\nContradiction found for neighborhood(s): {list(contradictory_neighborhoods.keys())}")
    print("This means no single elementary CA rule can generate the exact pattern.")
    print("We will generate candidate rules by resolving the contradiction for each possibility.")
    
    possible_rules = []
    
    # We found that (1,0,1) is the problematic neighborhood
    # First candidate: assume (1,0,1) -> 0
    rule_bits_1 = []
    rule_1_map = {}
    print("\nCase 1: Assume (1,0,1) -> 0")
    for i in range(7, -1, -1):
        neighborhood = tuple(map(int, f'{i:03b}'))
        if neighborhood == (1,0,1):
            output = 0
        else:
            # For non-contradictory ones, take the single observed value
            # Assume 0 for any that were truly never observed
            output = list(transitions.get(neighborhood, {0}))[0]
        rule_bits_1.append(str(output))
        rule_1_map[neighborhood] = output

    rule_1_str = "".join(rule_bits_1)
    rule_1_dec = int(rule_1_str, 2)
    possible_rules.append(rule_1_dec)
    print(f"Resulting rule map: {rule_1_map}")
    print(f"Rule binary: {rule_1_str}, Decimal: {rule_1_dec}")

    # Second candidate: assume (1,0,1) -> 1
    rule_bits_2 = []
    rule_2_map = {}
    print("\nCase 2: Assume (1,0,1) -> 1")
    for i in range(7, -1, -1):
        neighborhood = tuple(map(int, f'{i:03b}'))
        if neighborhood == (1,0,1):
            output = 1
        else:
            output = list(transitions.get(neighborhood, {0}))[0]
        rule_bits_2.append(str(output))
        rule_2_map[neighborhood] = output
        
    rule_2_str = "".join(rule_bits_2)
    rule_2_dec = int(rule_2_str, 2)
    possible_rules.append(rule_2_dec)
    print(f"Resulting rule map: {rule_2_map}")
    print(f"Rule binary: {rule_2_str}, Decimal: {rule_2_dec}")

    # Step 5: Sort and report the rules
    possible_rules.sort()
    
    final_answer_str = ",".join(map(str, possible_rules))
    print(f"\nThe possible rules, sorted in increasing order, are: {final_answer_str}")
    return possible_rules

solve_cellular_automaton()
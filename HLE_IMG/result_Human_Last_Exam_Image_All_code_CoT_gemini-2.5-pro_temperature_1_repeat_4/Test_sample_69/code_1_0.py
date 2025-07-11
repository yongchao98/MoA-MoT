import itertools

def solve_cellular_automaton():
    """
    Determines the possible elementary cellular automaton rules that could
    have generated the given visual pattern.
    """
    # 1. Transcribe the pattern from the image into a grid of 0s and 1s.
    grid = [
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0],
        [0,0,0,0,1,0,1,0,1,0,1,0,1,0,0,0,0],
        [0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0],
        [0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0],
        [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],
        [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]
    ]

    # 2. Extract transition constraints from the pattern.
    constraints = {}
    height = len(grid)
    width = len(grid[0])
    for t in range(height - 1):
        for i in range(1, width - 1):
            neighborhood = (grid[t][i-1], grid[t][i], grid[t][i+1])
            result = grid[t+1][i]
            if neighborhood in constraints and constraints[neighborhood] != result:
                # This should not happen for a valid ECA pattern.
                raise ValueError(f"Contradiction found for neighborhood {neighborhood}")
            constraints[neighborhood] = result

    # 3. Identify which parts of the rule are constrained and which are not.
    # The standard order of neighborhoods for defining the rule number.
    neighborhood_order = [
        (1,1,1), (1,1,0), (1,0,1), (1,0,0),
        (0,1,1), (0,1,0), (0,0,1), (0,0,0)
    ]

    rule_template = []
    unconstrained_indices = []
    for i, n_hood in enumerate(neighborhood_order):
        if n_hood in constraints:
            rule_template.append(constraints[n_hood])
        else:
            # Mark this position as a "don't care"
            rule_template.append(None)
            unconstrained_indices.append(i)
            
    # 4. Generate all possible rules by filling in the "don't care" positions.
    possible_rules = []
    num_unconstrained = len(unconstrained_indices)

    # Iterate through all combinations (e.g., (0,0,0), (0,0,1), ...) for the unconstrained bits.
    for fillings in itertools.product([0, 1], repeat=num_unconstrained):
        current_rule_bits = list(rule_template)
        for i, val in enumerate(fillings):
            # Fill in the "don't care" positions
            current_rule_bits[unconstrained_indices[i]] = val
        
        # Convert the list of bits to a binary string
        rule_binary_string = "".join(map(str, current_rule_bits))
        # Convert the binary string to its decimal rule number
        rule_decimal = int(rule_binary_string, 2)
        possible_rules.append(rule_decimal)

    # 5. Sort the rules and print the result.
    possible_rules.sort()
    print("The possible rules are:")
    print(','.join(map(str, possible_rules)))

solve_cellular_automaton()
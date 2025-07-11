import itertools

def find_possible_rules():
    """
    Finds possible elementary cellular automaton rules based on the first few rows of the given pattern.
    """
    # Grid transcription of the first 3 rows from the image (0=white, 1=black).
    # The pattern starts with a single black cell at the top.
    # The grid is wide enough to avoid boundary effects for the first few steps.
    grid_data = [
        "000000010000000",
        "000000111000000",
        "000001101100000",
    ]
    grid = [[int(c) for c in row] for row in grid_data]
    height = len(grid)
    width = len(grid[0])

    observed_constraints = {}

    # Extract constraints from the transitions between rows
    for t in range(height - 1):
        # We assume the grid is surrounded by 0s (quiescent background)
        current_row = [0] + grid[t] + [0]
        next_row = grid[t+1]
        
        for i in range(width):
            # The neighborhood consists of cells (i-1, i, i+1) in the original grid
            # which corresponds to (i, i+1, i+2) in the padded 'current_row'.
            neighborhood = tuple(current_row[i:i+3])
            result = next_row[i]

            if neighborhood in observed_constraints and observed_constraints[neighborhood] != result:
                # This should not happen in the first few clean rows.
                print(f"Error: Contradiction found for neighborhood {neighborhood}.")
                return

            observed_constraints[neighborhood] = result

    # All possible 3-cell neighborhoods in Wolfram order (from 111 down to 000)
    all_neighborhoods = list(itertools.product([1, 0], repeat=3))

    possible_rules = []

    # Iterate through all 256 possible rules
    for rule_num in range(256):
        is_consistent = True
        # Get the 8-bit binary representation of the rule number
        rule_bits = format(rule_num, '08b')
        
        # Create a mapping from neighborhood to outcome for the current rule
        rule_map = {all_neighborhoods[i]: int(rule_bits[i]) for i in range(8)}

        # Check if this rule is consistent with all observed constraints
        for neighborhood, result in observed_constraints.items():
            if rule_map[neighborhood] != result:
                is_consistent = False
                break
        
        if is_consistent:
            possible_rules.append(rule_num)
            
    # Print the final result
    print(','.join(map(str, sorted(possible_rules))))

find_possible_rules()
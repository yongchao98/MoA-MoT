import collections

def solve_cellular_automaton():
    """
    Identifies the elementary cellular automaton rules that could produce the given pattern.
    """
    # Step 1: Define the observed grid pattern from the image.
    # 0 represents a white square, 1 represents a black square.
    observed_grid = [
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0],
        [0,0,0,0,1,0,1,0,1,0,1,0,1,0,0,0,0],
        [0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0],
        [0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0],
        [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],
        [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1],
    ]

    # Step 2: Determine the rule constraints from the observed evolution.
    # The constraints dictionary will store {neighborhood_tuple: output_state}.
    constraints = {}
    height = len(observed_grid)
    width = len(observed_grid[0])

    for t in range(height - 1):
        current_row = observed_grid[t]
        next_row = observed_grid[t+1]
        # Iterate through each cell to find the neighborhood and its outcome.
        # Boundary is assumed to be all zeros.
        for i in range(1, width - 1):
            neighborhood = (current_row[i-1], current_row[i], current_row[i+1])
            output = next_row[i]

            if neighborhood in constraints:
                # If rule is already deduced, check for consistency.
                if constraints[neighborhood] != output:
                    # This case should not be reached for a valid CA pattern.
                    print(f"Error: Inconsistent rule for neighborhood {neighborhood}.")
                    return
            else:
                # Add the newly found transition rule to our constraints.
                constraints[neighborhood] = output

    # Step 3: Find all rules (0-255) that satisfy the discovered constraints.
    matching_rules = []
    for rule_number in range(256):
        is_a_match = True
        for neighborhood, required_output in constraints.items():
            # The standard ECA convention maps neighborhoods ('111' to '000')
            # to the bits of the rule number (bit 7 to bit 0).
            l, c, r = neighborhood
            bit_position = l * 4 + c * 2 + r
            
            # Check if the rule's output for this neighborhood matches the requirement.
            # (rule_number >> bit_position) & 1 extracts the relevant bit.
            if ((rule_number >> bit_position) & 1) != required_output:
                is_a_match = False
                break  # This rule is invalid, move to the next.
        
        if is_a_match:
            matching_rules.append(rule_number)

    # Step 4: Print the final sorted list of matching rule numbers.
    print(','.join(map(str, matching_rules)))

solve_cellular_automaton()
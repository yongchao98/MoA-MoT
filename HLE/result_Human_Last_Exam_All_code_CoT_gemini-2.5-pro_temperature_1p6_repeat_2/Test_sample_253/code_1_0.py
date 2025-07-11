def solve_rule_110():
    """
    This function simulates the Rule 110 cellular automaton for 20 generations
    starting from a single cell in state 1 and prints the final pattern.
    """
    # 1. Define simulation parameters
    generations = 20
    # The width must be large enough to contain the pattern's growth.
    # The pattern grows by at most 1 cell on each side per generation.
    width = 2 * generations + 1

    # 2. Define Rule 110. The key is the 3-cell neighborhood (left, center, right)
    # and the value is the next state of the center cell.
    rule110 = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0,
    }

    # 3. Set the initial state (generation 0): a single 1 in the center.
    current_gen = [0] * width
    current_gen[width // 2] = 1

    # 4. Run the simulation for the specified number of generations.
    for _ in range(generations):
        next_gen = [0] * width
        # Iterate through each cell to calculate its state in the next generation.
        # We can skip the very edges because their neighborhood will always be (0,0,0),
        # resulting in 0, which is already the default value in next_gen.
        for i in range(1, width - 1):
            # Get the 3-cell neighborhood as a tuple
            neighborhood = tuple(current_gen[i-1 : i+2])
            # Apply Rule 110 to determine the cell's next state
            next_gen[i] = rule110[neighborhood]
        
        # The newly calculated generation becomes the current one for the next iteration.
        current_gen = next_gen

    # 5. Print the final result.
    final_pattern_list = [str(cell) for cell in current_gen]
    final_pattern_str = "".join(final_pattern_list)
    print("The final pattern after 20 generations is:")
    print(final_pattern_str)

solve_rule_110()
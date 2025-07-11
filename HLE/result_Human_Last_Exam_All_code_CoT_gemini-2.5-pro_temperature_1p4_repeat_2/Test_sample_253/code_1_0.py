def solve_rule_110():
    """
    This function simulates the cellular automaton Rule 110 for 20 generations,
    starting from a single cell in state 1, and prints the final pattern.
    """
    # Parameters for the simulation
    num_generations = 20
    # The width needs to be large enough to contain the pattern.
    # The pattern grows by at most 1 cell on each side per generation.
    # Min width = 1 (initial) + 2 * num_generations = 41. We use a larger width for safety.
    width = 2 * num_generations + 21

    # The ruleset for Rule 110. The key is the 3-cell pattern (left, center, right),
    # and the value is the state of the center cell in the next generation.
    # The rule number 110 in binary is 01101110.
    rules = {
        '111': 0, '110': 1, '101': 1, '100': 0,
        '011': 1, '010': 1, '001': 1, '000': 0
    }

    # Initialize the grid: a line of cells, all in state 0,
    # except for a single cell in state 1 in the middle.
    line = [0] * width
    line[width // 2] = 1

    # Run the simulation for the specified number of generations
    for _ in range(num_generations):
        # Create a new line for the next generation, initialized to all zeros.
        next_line = [0] * width
        
        # Iterate through each cell to calculate its next state.
        # We skip the very edges as they are assumed to be in an infinite field of 0s.
        for i in range(1, width - 1):
            # Get the 3-cell neighborhood pattern as a string
            pattern = f"{line[i-1]}{line[i]}{line[i+1]}"
            
            # Look up the new state from the rule dictionary and update the next line
            next_line[i] = rules[pattern]
            
        # The new line becomes the current line for the next iteration
        line = next_line

    # Prepare the final output
    # Convert the list of integers into a single string
    result_string = "".join(map(str, line))
    
    # Trim leading and trailing zeros to get only the active pattern
    trimmed_result_string = result_string.strip('0')

    # Format the result by joining each digit with a space.
    # The problem asks to "output each number in the final equation" (likely meaning generation),
    # so we will print each digit separated by a space.
    if not trimmed_result_string:
         # Handles the case where the result is empty or all zeros
        print("0")
    else:
        final_output = " ".join(list(trimmed_result_string))
        print(final_output)

solve_rule_110()
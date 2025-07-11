import numpy as np

def solve_grid():
    """
    Solves the grid transformation puzzle by applying a derived rule.
    The rule determines the output of a cell based on its input value (C)
    and the sum of its 8 neighbors (S8).
    """
    
    # Test Input Grid 6
    input_grid = np.array([
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 0, 1],
        [1, 0, 1, 1, 0]
    ])
    
    # The "equation" or rule discovered from the examples.
    # It maps a tuple (C, S8) to the output value O.
    # C = cell's original value
    # S8 = sum of the 8 neighbors
    rule_map = {
        # Rules for when the cell value (C) is 0
        (0, 0): 0, (0, 1): 1, (0, 2): 1, (0, 3): 1, (0, 4): 0, (0, 5): 1, (0, 8): 0,
        # Rules for when the cell value (C) is 1
        (1, 1): 0, (1, 2): 1, (1, 3): 0, (1, 4): 0, (1, 5): 1, (1, 6): 0
    }

    # Get grid dimensions
    rows, cols = input_grid.shape
    
    # Initialize the output grid
    output_grid = np.zeros_like(input_grid)

    # Apply the rule to each cell of the input grid
    for r in range(rows):
        for c in range(cols):
            # Get the cell's current value (C)
            c_val = input_grid[r, c]
            
            # Calculate the sum of its 8 neighbors (S8)
            # Define slices for the 3x3 neighborhood
            r_slice = slice(max(0, r - 1), min(rows, r + 2))
            c_slice = slice(max(0, c - 1), min(cols, c + 2))
            # Sum the neighborhood and subtract the central cell's value
            s8_val = np.sum(input_grid[r_slice, c_slice]) - c_val
            
            # Look up the output value from the rule map
            output_val = rule_map.get((c_val, s8_val), 0) # Default to 0 if rule not found
            
            output_grid[r, c] = output_val

    # Print the final grid as a 25-digit binary string
    final_string = "".join(map(str, output_grid.flatten()))
    print(final_string)

solve_grid()
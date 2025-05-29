import numpy as np

def transform_grid(input_grid):
    # Find the positions of non-zero numbers
    positions = np.argwhere(input_grid > 0)
    output_blocks = []
    
    for pos in positions:
        value = input_grid[pos[0], pos[1]]
        # Create a 3x3 block for each non-zero number
        block = np.array([[0, value, 0],
                          [value, value, value],
                          [0, value, 0]])
        output_blocks.append(block)
    
    # Combine the blocks into the final output grid
    output_grid = np.vstack(output_blocks)
    return output_grid

# Test input grid
input_grid = np.array([
    [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0, 5, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
])

output_grid = transform_grid(input_grid)
print(output_grid)
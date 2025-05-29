import numpy as np

def transform_grid(input_grid):
    # Convert input to numpy array
    grid = np.array([[int(x) for x in row] for row in input_grid])
    
    # Create 8x8 output grid
    output = np.zeros((8, 8), dtype=int)
    
    # Top-left quadrant: Flip left-right
    top_left = np.fliplr(grid)
    output[0:4, 0:4] = top_left
    
    # Top-right quadrant: Original grid
    output[0:4, 4:8] = grid
    
    # Bottom half: Mirror of top half with vertical flip
    output[4:8, :] = np.flipud(output[0:4, :])
    
    return output

# Test input
input_grid = [
    [4, 6, 4, 2],
    [2, 5, 2, 4],
    [6, 5, 6, 6],
    [6, 5, 6, 6]
]

result = transform_grid(input_grid)

# Print the result in the required format
for row in result:
    print(' '.join(map(str, row)))
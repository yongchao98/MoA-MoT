import numpy as np

def extract_pattern(input_grid):
    # Convert input grid to numpy array for easier manipulation
    grid = np.array(input_grid)
    
    # Extract the first two rows and first two columns
    output_grid = grid[:2, :2]
    
    return output_grid.tolist()

# Test input grid
input_grid = [
    [3, 3, 8],
    [3, 3, 3],
    [8, 3, 8]
]

# Find the pattern and print the output grid
output_grid = extract_pattern(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))
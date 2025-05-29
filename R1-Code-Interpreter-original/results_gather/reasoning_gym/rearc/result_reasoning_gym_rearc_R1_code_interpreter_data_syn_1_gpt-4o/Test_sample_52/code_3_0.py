def extract_diagonal_pattern(input_grid):
    # Convert input grid to numpy array for easier manipulation
    grid = np.array(input_grid)
    
    # Extract the diagonal elements
    diagonal = np.diagonal(grid)
    
    # Create the output grid based on the diagonal elements
    output_grid = [[element] for element in diagonal]
    
    return output_grid

# Test input grid
input_grid = [
    [3, 3, 8],
    [3, 3, 3],
    [8, 3, 8]
]

# Find the pattern and print the output grid
output_grid = extract_diagonal_pattern(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))
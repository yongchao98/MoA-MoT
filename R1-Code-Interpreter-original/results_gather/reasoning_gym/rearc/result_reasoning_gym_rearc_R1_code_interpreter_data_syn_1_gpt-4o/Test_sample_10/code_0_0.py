def extract_subgrid(input_grid):
    # Extract the non-border elements
    return [row[1:-1] for row in input_grid[1:-1]]

# Test input grid
input_grid = [
    [4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4],
    [4, 6, 8, 6, 4],
    [4, 4, 4, 4, 4],
    [4, 4, 4, 4, 4]
]

# Extract the subgrid
output_grid = extract_subgrid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
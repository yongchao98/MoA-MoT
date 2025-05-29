def extract_subgrid(input_grid):
    # Find the bounds of the non-zero subgrid
    min_row, max_row = len(input_grid), 0
    min_col, max_col = len(input_grid[0]), 0
    
    for r in range(len(input_grid)):
        for c in range(len(input_grid[0])):
            if input_grid[r][c] != 0:
                min_row = min(min_row, r)
                max_row = max(max_row, r)
                min_col = min(min_col, c)
                max_col = max(max_col, c)
    
    # Extract the subgrid
    subgrid = []
    for r in range(min_row, max_row + 1):
        subgrid.append(input_grid[r][min_col:max_col + 1])
    
    # Transform the subgrid if needed (e.g., add counts or sums)
    # For now, we assume no additional transformation is needed
    return subgrid

# Test input grid
test_input = [
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 5, 2, 0, 0],
    [0, 0, 2, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
]

# Extract and print the output grid
output_grid = extract_subgrid(test_input)
print(output_grid)
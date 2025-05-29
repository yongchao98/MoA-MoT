def extract_subgrid(input_grid):
    # Find the bounds of the non-5 numbers
    min_row, max_row = len(input_grid), -1
    min_col, max_col = len(input_grid[0]), -1
    
    for r in range(len(input_grid)):
        for c in range(len(input_grid[0])):
            if input_grid[r][c] != 5:
                min_row = min(min_row, r)
                max_row = max(max_row, r)
                min_col = min(min_col, c)
                max_col = max(max_col, c)
    
    # Calculate the starting indices for the 4x4 subgrid
    start_row = max_row - 3
    start_col = max_col - 3
    
    # Ensure the subgrid is within bounds
    if start_row < min_row:
        start_row = min_row
    if start_col < min_col:
        start_col = min_col
    
    # Extract the 4x4 subgrid
    output_grid = []
    for r in range(start_row, start_row + 4):
        output_grid.append(input_grid[r][start_col:start_col + 4])
    
    return output_grid

# Test input grid
input_grid = [
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 2, 2, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 0, 5, 5, 5, 5, 5, 1, 1, 5, 5],
    [5, 5, 0, 0, 5, 5, 5, 5, 5, 1, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
]

output_grid = extract_subgrid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
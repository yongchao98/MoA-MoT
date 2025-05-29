def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Get the dimensions of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Iterate over the grid to find regions and apply the rule
    for r in range(rows):
        for c in range(cols):
            # Check if the current cell is part of a region
            if input_grid[r][c] != 0:
                # Check the boundaries of the region
                if (r == 0 or input_grid[r-1][c] == 0) and (r < rows - 1 and input_grid[r+1][c] == input_grid[r][c]):
                    # Top boundary
                    output_grid[r][c] = 2
                if (r == rows - 1 or input_grid[r+1][c] == 0) and (r > 0 and input_grid[r-1][c] == input_grid[r][c]):
                    # Bottom boundary
                    output_grid[r][c] = 2
                if (c == 0 or input_grid[r][c-1] == 0) and (c < cols - 1 and input_grid[r][c+1] == input_grid[r][c]):
                    # Left boundary
                    output_grid[r][c] = 2
                if (c == cols - 1 or input_grid[r][c+1] == 0) and (c > 0 and input_grid[r][c-1] == input_grid[r][c]):
                    # Right boundary
                    output_grid[r][c] = 2
    
    return output_grid

# Test input grid
input_grid = [
    [5, 0, 0, 0, 0, 0, 5, 5],
    [5, 0, 5, 5, 5, 0, 5, 5],
    [5, 0, 0, 0, 0, 0, 5, 5],
    [5, 0, 5, 5, 5, 0, 5, 5],
    [5, 0, 0, 0, 0, 0, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
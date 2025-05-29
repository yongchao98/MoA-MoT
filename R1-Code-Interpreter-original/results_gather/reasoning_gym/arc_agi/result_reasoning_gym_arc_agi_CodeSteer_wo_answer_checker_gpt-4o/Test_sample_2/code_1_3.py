def transform_grid(input_grid):
    # Convert the input grid to a list of lists for easier manipulation
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    rows = len(grid)
    cols = len(grid[0])
    
    # Create a copy of the grid to store the output
    output_grid = [row[:] for row in grid]
    
    # Iterate over the grid, skipping the border to avoid index errors
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            # Check if the current position is the center of a 3x3 block of '5's
            if (grid[i][j] == 5 and
                grid[i-1][j-1] == 5 and grid[i-1][j] == 5 and grid[i-1][j+1] == 5 and
                grid[i][j-1] == 5 and grid[i][j+1] == 5 and
                grid[i+1][j-1] == 5 and grid[i+1][j] == 5 and grid[i+1][j+1] == 5):
                # Replace the 3x3 block with '0's in the output grid
                output_grid[i-1][j-1] = 0
                output_grid[i-1][j] = 0
                output_grid[i-1][j+1] = 0
                output_grid[i][j-1] = 0
                output_grid[i][j] = 0
                output_grid[i][j+1] = 0
                output_grid[i+1][j-1] = 0
                output_grid[i+1][j] = 0
                output_grid[i+1][j+1] = 0
    
    # Convert the output grid back to the required format
    output_grid_str = '\n'.join(' '.join(map(str, row)) for row in output_grid)
    return output_grid_str

# Test input grid
test_input = """
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 5 5 4 5 5 5 4 4 4 5 5 5 4 4 5 5 5 4
4 4 5 4 5 4 5 4 5 4 4 4 5 4 5 4 4 5 4 5 4
4 4 5 5 5 4 5 5 5 4 4 4 5 5 5 4 4 5 5 5 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 4 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 4 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 4 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 4 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
"""

# Get the output grid
output = transform_grid(test_input)
print(f"<<<{output}>>>")
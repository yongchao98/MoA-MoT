def transform_grid(input_grid):
    # Convert the input grid into a list of lists for easier manipulation
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    
    # Create a copy of the grid to store the output
    output_grid = [row[:] for row in grid]
    
    # Iterate over the grid to find 3x3 blocks of '5's
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            # Check if the current position is the center of a 3x3 block of '5's
            if (grid[i-1][j-1] == grid[i-1][j] == grid[i-1][j+1] == 
                grid[i][j-1] == grid[i][j] == grid[i][j+1] == 
                grid[i+1][j-1] == grid[i+1][j] == grid[i+1][j+1] == 5):
                
                # Change the outer '5's to '0'
                output_grid[i-1][j-1] = 0
                output_grid[i-1][j] = 0
                output_grid[i-1][j+1] = 0
                output_grid[i][j-1] = 0
                # output_grid[i][j] remains 5
                output_grid[i][j+1] = 0
                output_grid[i+1][j-1] = 0
                output_grid[i+1][j] = 0
                output_grid[i+1][j+1] = 0
    
    # Convert the output grid back to the required format
    return '\n'.join(' '.join(map(str, row)) for row in output_grid)

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

# Transform the test input grid
output_grid = transform_grid(test_input)
print(output_grid)
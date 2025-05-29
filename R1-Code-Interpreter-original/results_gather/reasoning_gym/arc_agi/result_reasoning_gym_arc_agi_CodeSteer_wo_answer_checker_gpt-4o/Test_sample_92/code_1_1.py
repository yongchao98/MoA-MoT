def transform_grid(input_grid):
    # Convert the input grid into a 2D list
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    
    # Create a copy of the grid for the output
    output_grid = [row[:] for row in grid]
    
    # Function to extend a number in a specific pattern
    def extend_pattern(r, c, num):
        # Example pattern: extend in a cross shape
        for i in range(max(0, r-1), min(rows, r+2)):
            for j in range(max(0, c-1), min(cols, c+2)):
                if grid[i][j] == 3:  # Only replace default values
                    output_grid[i][j] = num
    
    # Iterate over each cell in the grid
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] != 3:  # If the cell is not the default value
                num = grid[r][c]
                extend_pattern(r, c, num)
    
    # Convert the output grid back to the required format
    output_str = '\n'.join(' '.join(map(str, row)) for row in output_grid)
    return output_str

# Test input grid
test_input = """
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 7 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 7 7 7 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 9 7 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 6 1 1 1 3 3 3 3 3 3 5 9 3 3 3 3 3 3 3
3 3 3 3 3 7 5 3 3 3 3 3 3 3 3 7 6 6 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 6 3 6 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 6 6 6 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 6 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
"""

# Get the output grid
output_grid = transform_grid(test_input)
print("<<<" + output_grid + ">>>")
def transform_grid(input_grid):
    # Convert the input grid into a 2D list
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    rows, cols = len(grid), len(grid[0])
    
    # Find the largest number in the grid
    largest_number = max(max(row) for row in grid)
    
    # Determine the new number to be used for replacement
    new_number = largest_number + 1
    
    # Create a copy of the grid to modify
    output_grid = [row[:] for row in grid]
    
    # Calculate the middle row and column
    mid_row = rows // 2
    mid_col = cols // 2
    
    # Apply the transformation rule: replace the middle row and column with the new number
    for c in range(cols):
        output_grid[mid_row][c] = new_number
    for r in range(rows):
        output_grid[r][mid_col] = new_number
    
    # Convert the output grid back to the required format
    output_grid_str = '\n'.join(' '.join(map(str, row)) for row in output_grid)
    return output_grid_str

# Test input grid
test_input = """
2 2 2 3 5 2 5 3 5 5 5 3 2 2 2 3 5 2 2 3 2 5 2
2 2 2 3 5 5 5 3 5 5 5 3 2 5 5 3 5 2 5 3 5 5 2
2 5 5 3 2 2 5 3 5 2 5 3 2 5 2 3 5 5 2 3 2 5 2
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
5 5 2 3 2 2 2 3 5 2 2 3 5 2 5 3 2 2 2 3 2 5 2
2 2 5 3 2 5 2 3 5 2 2 3 5 2 2 3 5 5 2 3 5 2 2
2 2 2 3 5 5 2 3 5 2 2 3 2 2 2 3 5 2 2 3 2 5 2
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
2 5 5 3 2 5 2 3 5 2 5 3 2 2 5 3 5 2 2 3 5 2 2
5 2 2 3 2 2 2 3 2 2 5 3 5 5 5 3 5 5 2 3 2 2 2
5 2 5 3 2 2 2 3 5 2 5 3 5 5 2 3 2 2 2 3 2 2 5
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
2 2 2 3 5 2 2 3 5 5 2 3 5 2 5 3 5 2 2 3 5 2 2
5 2 2 3 5 5 5 3 2 2 2 3 5 5 2 3 5 2 5 3 5 2 2
5 2 2 3 2 2 5 3 5 5 2 3 5 5 2 3 2 2 5 3 5 5 2
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
5 2 2 3 5 2 5 3 5 2 2 3 2 5 2 3 1 1 2 2 2 2 3
2 5 2 3 5 2 2 3 2 5 5 3 2 5 5 3 1 1 2 2 2 2 3
2 2 2 3 2 2 2 3 5 5 5 3 2 2 2 3 1 1 1 1 1 1 3
3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 1 1 1 1 3
2 5 2 3 2 2 2 3 2 2 2 3 5 5 2 3 2 2 2 2 1 1 3
5 2 2 3 2 2 2 3 2 2 5 3 5 5 5 3 2 2 2 2 1 1 3
5 2 5 3 5 2 2 3 2 5 5 3 2 5 2 3 3 3 3 3 3 3 3
"""

# Get the output grid
output_grid = transform_grid(test_input)
print(output_grid)
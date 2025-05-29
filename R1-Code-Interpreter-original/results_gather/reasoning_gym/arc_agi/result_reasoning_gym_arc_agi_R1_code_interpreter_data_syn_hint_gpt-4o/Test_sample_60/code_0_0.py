def transform_grid(input_grid):
    # Convert the input grid into a 2D list
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    
    # Find the largest number in the grid, excluding 3
    largest_number = max(num for row in grid for num in row if num != 3)
    
    # Create a copy of the grid to modify
    output_grid = [row[:] for row in grid]
    
    # Determine the middle section to replace with the largest number
    # Assuming the middle section is a 3x3 block in the center of the grid
    rows, cols = len(grid), len(grid[0])
    mid_row_start, mid_row_end = rows // 2 - 1, rows // 2 + 2
    mid_col_start, mid_col_end = cols // 2 - 1, cols // 2 + 2
    
    # Replace the middle section with the largest number
    for i in range(mid_row_start, mid_row_end):
        for j in range(mid_col_start, mid_col_end):
            output_grid[i][j] = largest_number
    
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
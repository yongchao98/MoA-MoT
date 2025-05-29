def transform_grid(input_grid):
    # Convert the input grid to a list of lists for easier manipulation
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]
    
    # Find the distinct region (the region surrounded by a different number)
    distinct_value = None
    surrounding_value = None
    
    # Determine the distinct value and surrounding value
    for row in grid:
        for value in row:
            if distinct_value is None:
                distinct_value = value
            elif value != distinct_value:
                surrounding_value = value
                break
        if surrounding_value is not None:
            break
    
    # Find the bounds of the distinct region
    top, bottom, left, right = None, None, None, None
    for i, row in enumerate(grid):
        for j, value in enumerate(row):
            if value == distinct_value:
                if top is None:
                    top = i
                bottom = i
                if left is None or j < left:
                    left = j
                if right is None or j > right:
                    right = j
    
    # Insert a line of '1's in the second row of the distinct region
    for j in range(left, right + 1):
        grid[top + 1][j] = 1
    
    # Convert the grid back to the required output format
    output_grid = '\n'.join(' '.join(map(str, row)) for row in grid)
    return output_grid

# Test input grid
test_input = """
6 6 6 6 6 6 6 6 6 6 6
6 6 6 6 6 6 6 6 6 6 6
0 0 0 0 0 0 0 0 6 6 6
0 6 6 6 0 6 6 0 6 6 6
0 0 0 0 0 0 0 0 6 6 6
6 6 6 6 6 6 6 6 6 6 6
6 6 6 6 6 6 6 6 6 6 6
6 6 6 6 6 6 6 6 6 6 6
0 0 0 0 0 0 6 6 6 6 6
0 6 6 6 6 0 6 6 6 6 6
6 6 0 6 0 0 6 6 6 6 6
"""

# Transform the test input grid
output_grid = transform_grid(test_input)
print(output_grid)
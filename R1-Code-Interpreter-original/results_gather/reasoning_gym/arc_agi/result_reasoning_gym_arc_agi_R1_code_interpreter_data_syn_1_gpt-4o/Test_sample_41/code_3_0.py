def transform_grid(input_grid):
    # Convert the input grid to a list of lists for easier manipulation
    grid = [list(map(int, row.split())) for row in input_grid.strip().split('\n')]

    # Apply the transformation rule
    # Ensure that rows 8 and 9 are identical
    grid[5] = grid[6] = grid[7]
    grid[8] = grid[9] = grid[10]

    # Convert the grid back to the required output format
    output_grid = '\n'.join(' '.join(map(str, row)) for row in grid)
    return output_grid

# Test input grid
input_grid = """
4 4 7 5 2 7 9 8 8 9 7 2 5 7 4 4
4 5 5 7 7 7 8 8 8 8 7 7 7 5 5 4
7 5 8 2 9 8 8 3 3 3 8 9 2 8 5 7
5 7 2 5 8 8 5 3 3 3 8 8 5 2 7 5
2 7 3 3 3 3 8 3 3 3 4 7 8 9 7 2
7 7 3 3 3 3 8 0 0 8 4 4 8 8 7 7
9 8 3 3 3 3 0 0 0 0 8 8 5 8 8 9
8 8 3 3 3 3 0 4 4 0 0 8 5 5 8 8
8 8 5 5 8 0 0 4 4 0 0 8 5 5 8 8
9 8 8 5 8 8 0 0 0 0 8 8 5 8 8 9
7 7 8 8 4 4 8 0 0 8 4 4 8 8 7 7
2 7 9 8 7 4 8 8 8 8 4 7 8 9 7 2
5 7 2 5 8 8 5 5 5 5 8 8 5 2 7 5
7 5 8 2 9 8 8 5 5 8 8 9 2 8 5 7
4 5 5 7 7 7 8 8 8 8 7 7 7 5 5 4
4 4 7 5 2 7 9 8 8 9 7 2 5 7 4 4
"""

# Apply the transformation and print the output grid
output_grid = transform_grid(input_grid)
print(output_grid)
def remove_border(grid):
    # Check if the grid is too small to have a border
    if len(grid) <= 2 or len(grid[0]) <= 2:
        return []  # Return an empty grid if there's no inner grid

    # Remove the first and last row
    inner_grid = grid[1:-1]

    # Remove the first and last column from each row
    inner_grid = [row[1:-1] for row in inner_grid]

    return inner_grid

# Test the function with the provided input grid
input_grid = [
    [9, 9, 9, 9, 9, 9],
    [9, 7, 5, 7, 5, 7],
    [9, 5, 7, 5, 7, 5],
    [9, 7, 5, 7, 5, 7],
    [9, 5, 7, 5, 7, 5],
    [9, 7, 5, 7, 5, 7]
]

output_grid = remove_border(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
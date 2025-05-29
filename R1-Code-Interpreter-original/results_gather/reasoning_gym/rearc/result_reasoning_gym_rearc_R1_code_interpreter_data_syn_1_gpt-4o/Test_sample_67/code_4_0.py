def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid to output grid

    # A helper function to perform flood fill and mark the surrounded zeros
    def flood_fill(i, j):
        if i < 0 or i >= rows or j < 0 or j >= cols:
            return False
        if input_grid[i][j] != 0:
            return True
        if output_grid[i][j] == 1:
            return True

        # Temporarily mark the cell to avoid revisiting
        output_grid[i][j] = 1

        # Check all four directions
        up = flood_fill(i - 1, j)
        down = flood_fill(i + 1, j)
        left = flood_fill(i, j - 1)
        right = flood_fill(i, j + 1)

        # If any direction is not surrounded, this cell is not surrounded
        if not (up and down and left and right):
            output_grid[i][j] = 0
            return False

        return True

    # Apply flood fill to all zeros
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 0 and output_grid[i][j] != 1:
                flood_fill(i, j)

    return output_grid

# Test input grid
input_grid = [
    [3, 3, 0, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 0],
    [3, 3, 3, 3, 0, 0, 0],
    [3, 0, 0, 3, 0, 0, 0],
    [0, 3, 0, 3, 0, 0, 0],
    [0, 3, 0, 0, 0, 3, 3],
    [3, 3, 0, 3, 0, 3, 3]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    # Define a function to check and replace patterns
    def check_and_replace(r, c):
        # Check horizontal pattern
        if c <= cols - 3 and input_grid[r][c] == input_grid[r][c + 1] == input_grid[r][c + 2]:
            output_grid[r][c] = 5
            output_grid[r][c + 1] = 5
            output_grid[r][c + 2] = 5

        # Check vertical pattern
        if r <= rows - 3 and input_grid[r][c] == input_grid[r + 1][c] == input_grid[r + 2][c]:
            output_grid[r][c] = 5
            output_grid[r + 1][c] = 5
            output_grid[r + 2][c] = 5

        # Check diagonal pattern (top-left to bottom-right)
        if r <= rows - 3 and c <= cols - 3 and input_grid[r][c] == input_grid[r + 1][c + 1] == input_grid[r + 2][c + 2]:
            output_grid[r][c] = 5
            output_grid[r + 1][c + 1] = 5
            output_grid[r + 2][c + 2] = 5

        # Check diagonal pattern (top-right to bottom-left)
        if r <= rows - 3 and c >= 2 and input_grid[r][c] == input_grid[r + 1][c - 1] == input_grid[r + 2][c - 2]:
            output_grid[r][c] = 5
            output_grid[r + 1][c - 1] = 5
            output_grid[r + 2][c - 2] = 5

    # Apply the check and replace function to each cell
    for r in range(rows):
        for c in range(cols):
            check_and_replace(r, c)

    return output_grid

# Test input grid
input_grid = [
    [3, 9, 3, 3, 0, 3, 2, 4, 3, 3, 6, 3, 3, 3, 6],
    [3, 3, 3, 5, 3, 3, 7, 3, 3, 2, 6, 3, 3, 0, 0],
    [3, 6, 3, 3, 3, 6, 0, 3, 2, 7, 3, 2, 3, 3, 3],
    [7, 3, 7, 3, 0, 0, 3, 5, 0, 3, 3, 3, 3, 3, 3],
    [3, 4, 2, 9, 3, 0, 3, 3, 2, 7, 3, 5, 2, 3, 3],
    [3, 7, 3, 3, 5, 3, 3, 3, 0, 2, 4, 3, 3, 3, 3],
    [0, 3, 0, 3, 3, 3, 3, 4, 7, 3, 3, 2, 3, 0, 2],
    [3, 3, 9, 9, 3, 6, 7, 3, 0, 7, 3, 4, 0, 4, 3],
    [7, 7, 7, 0, 3, 2, 7, 3, 5, 3, 9, 7, 3, 5, 3],
    [3, 3, 3, 3, 9, 3, 6, 3, 0, 2, 3, 0, 4, 3, 3],
    [3, 7, 5, 3, 3, 3, 6, 3, 0, 3, 3, 7, 3, 4, 9],
    [3, 3, 3, 3, 3, 2, 2, 7, 3, 3, 4, 2, 3, 3, 9],
    [3, 3, 9, 3, 2, 3, 3, 3, 6, 2, 2, 3, 4, 3, 3],
    [3, 3, 7, 7, 0, 2, 3, 2, 3, 3, 6, 2, 3, 5, 3],
    [4, 2, 9, 2, 6, 7, 7, 3, 7, 6, 2, 7, 3, 4, 0]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
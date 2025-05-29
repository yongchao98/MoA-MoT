def expand_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    # Find the unique number that needs to be expanded
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] != input_grid[0][0]:  # Assuming the first element is the common one
                unique_number = input_grid[i][j]
                # Expand horizontally
                if j > 0:
                    output_grid[i][j-1] = unique_number
                if j < cols - 1:
                    output_grid[i][j+1] = unique_number
                # Expand vertically
                if i > 0:
                    output_grid[i-1][j] = unique_number
                if i < rows - 1:
                    output_grid[i+1][j] = unique_number

    return output_grid

# Test input grid
input_grid = [
    [1, 2, 2, 2, 2, 2],
    [1, 2, 2, 2, 2, 2],
    [1, 2, 2, 2, 2, 2],
    [1, 2, 2, 7, 2, 2],
    [1, 2, 2, 2, 2, 2],
    [1, 2, 2, 2, 2, 2],
    [1, 2, 2, 2, 2, 2]
]

output_grid = expand_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
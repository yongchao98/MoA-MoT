def transform_grid(input_grid):
    n = len(input_grid)
    m = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Create a copy of the input grid

    # Increment the diagonal elements
    for i in range(min(n, m)):
        output_grid[i][i] += 1

    # Set the last row to 1s, except the last element
    output_grid[-1] = [1] * (m - 1) + [input_grid[-1][-1]]

    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3],
    [8, 8, 8, 8, 8, 8, 8, 8, 8, 3]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
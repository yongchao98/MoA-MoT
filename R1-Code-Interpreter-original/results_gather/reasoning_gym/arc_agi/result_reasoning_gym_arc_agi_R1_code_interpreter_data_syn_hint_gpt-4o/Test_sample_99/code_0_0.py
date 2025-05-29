def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    # Identify the numbers for cross and diagonal patterns
    cross_number = 9
    diagonal_number = 5

    # Apply cross pattern
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == cross_number:
                # Create a cross pattern
                if i > 0:
                    output_grid[i-1][j] = cross_number
                if i < rows - 1:
                    output_grid[i+1][j] = cross_number
                if j > 0:
                    output_grid[i][j-1] = cross_number
                if j < cols - 1:
                    output_grid[i][j+1] = cross_number

    # Apply diagonal pattern
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == diagonal_number:
                # Create a diagonal pattern
                if i > 0 and j > 0:
                    output_grid[i-1][j-1] = diagonal_number
                if i > 0 and j < cols - 1:
                    output_grid[i-1][j+1] = diagonal_number
                if i < rows - 1 and j > 0:
                    output_grid[i+1][j-1] = diagonal_number
                if i < rows - 1 and j < cols - 1:
                    output_grid[i+1][j+1] = diagonal_number

    return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 9, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [2, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 9, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))
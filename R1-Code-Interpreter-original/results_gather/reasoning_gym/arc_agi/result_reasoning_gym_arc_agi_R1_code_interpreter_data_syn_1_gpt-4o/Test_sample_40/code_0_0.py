def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [row[:] for row in input_grid]  # Copy the input grid

    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 4:
                # Add 3s and 1s around the 4s
                if i > 0 and output_grid[i-1][j] == 0:
                    output_grid[i-1][j] = 3
                if i < rows - 1 and output_grid[i+1][j] == 0:
                    output_grid[i+1][j] = 3
                if j > 0 and output_grid[i][j-1] == 0:
                    output_grid[i][j-1] = 3
                if j < cols - 1 and output_grid[i][j+1] == 0:
                    output_grid[i][j+1] = 3
                # Diagonal positions for 1s
                if i > 0 and j > 0 and output_grid[i-1][j-1] == 0:
                    output_grid[i-1][j-1] = 1
                if i > 0 and j < cols - 1 and output_grid[i-1][j+1] == 0:
                    output_grid[i-1][j+1] = 1
                if i < rows - 1 and j > 0 and output_grid[i+1][j-1] == 0:
                    output_grid[i+1][j-1] = 1
                if i < rows - 1 and j < cols - 1 and output_grid[i+1][j+1] == 0:
                    output_grid[i+1][j+1] = 1

    return output_grid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0],
    [0, 0, 4, 4, 0, 0, 0, 0, 0, 0, 6, 4, 4, 0, 0],
    [0, 0, 4, 0, 4, 0, 0, 0, 0, 1, 4, 0, 4, 0, 0],
    [0, 0, 6, 4, 4, 0, 0, 0, 0, 3, 4, 4, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 4, 4, 6, 0, 0, 0, 0, 0, 4, 4, 0],
    [0, 0, 0, 0, 4, 0, 4, 0, 0, 0, 0, 4, 0, 4, 0],
    [0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 4, 4, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
]

output_grid = transform_grid(input_grid)

for row in output_grid:
    print(' '.join(map(str, row)))
def transform_grid(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output_grid = [[0] * cols for _ in range(rows)]

    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 8:
                # Check surrounding cells to decide if it should be 1
                if (i > 0 and input_grid[i-1][j] != 8) or (i < rows-1 and input_grid[i+1][j] != 8) or \
                   (j > 0 and input_grid[i][j-1] != 8) or (j < cols-1 and input_grid[i][j+1] != 8):
                    output_grid[i][j] = 1
                else:
                    output_grid[i][j] = 8
            else:
                output_grid[i][j] = input_grid[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 0, 8, 8, 4, 8, 8, 8, 9, 8, 2, 8, 8, 5],
    [8, 6, 5, 8, 0, 8, 4, 4, 8, 0, 9, 7, 8, 8, 8, 8],
    [7, 8, 8, 7, 3, 9, 1, 8, 4, 5, 7, 7, 8, 8, 8, 8],
    [0, 1, 8, 8, 8, 8, 8, 0, 8, 8, 8, 7, 4, 9, 3, 8],
    [8, 7, 8, 6, 9, 5, 0, 5, 8, 8, 0, 4, 8, 2, 9, 8],
    [3, 8, 8, 8, 4, 1, 8, 8, 8, 3, 8, 9, 0, 1, 6, 9],
    [2, 8, 2, 8, 8, 8, 8, 8, 0, 8, 8, 8, 8, 8, 0, 8],
    [3, 3, 1, 0, 8, 8, 8, 4, 8, 8, 3, 3, 3, 1, 8, 7],
    [8, 8, 8, 4, 8, 8, 8, 8, 8, 6, 8, 6, 3, 4, 8, 8],
    [2, 1, 8, 3, 2, 8, 8, 8, 8, 8, 8, 2, 8, 5, 3, 8],
    [7, 2, 8, 8, 8, 8, 9, 4, 1, 8, 8, 3, 8, 8, 8, 5],
    [8, 2, 4, 1, 7, 7, 8, 2, 9, 0, 5, 3, 7, 8, 4, 8],
    [8, 8, 3, 8, 5, 8, 9, 8, 1, 3, 8, 4, 8, 2, 8, 0],
    [6, 8, 8, 8, 2, 8, 8, 4, 8, 8, 0, 8, 0, 8, 9, 8],
    [0, 7, 8, 8, 8, 7, 8, 8, 6, 7, 1, 8, 5, 8, 2, 8],
    [1, 9, 1, 5, 9, 8, 8, 8, 1, 7, 8, 4, 8, 8, 9, 8]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))
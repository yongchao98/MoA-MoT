def should_replace_with_one(grid, i, j):
    # Check the surrounding elements of grid[i][j]
    rows, cols = len(grid), len(grid[0])
    neighbors = [
        (i-1, j-1), (i-1, j), (i-1, j+1),
        (i, j-1),           (i, j+1),
        (i+1, j-1), (i+1, j), (i+1, j+1)
    ]
    
    # Count how many neighbors are '8'
    eight_count = 0
    for ni, nj in neighbors:
        if 0 <= ni < rows and 0 <= nj < cols and grid[ni][nj] == 8:
            eight_count += 1
    
    # Example pattern: if '8' is surrounded by at least 5 other '8's, keep it as '8'
    if eight_count >= 5:
        return False
    return True

def transform_grid(input_grid):
    output_grid = []
    for i in range(len(input_grid)):
        output_row = []
        for j in range(len(input_grid[0])):
            if input_grid[i][j] == 8:
                if should_replace_with_one(input_grid, i, j):
                    output_row.append(1)
                else:
                    output_row.append(8)
            else:
                output_row.append(input_grid[i][j])
        output_grid.append(output_row)
    return output_grid

# Define the input grid
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

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
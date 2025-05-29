def extract_subgrid(grid):
    min_row, max_row = len(grid), -1
    min_col, max_col = len(grid[0]), -1

    # Find the bounds of the non-zero elements
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if grid[i][j] != 0:
                if i < min_row:
                    min_row = i
                if i > max_row:
                    max_row = i
                if j < min_col:
                    min_col = j
                if j > max_col:
                    max_col = j

    # Extract the subgrid
    subgrid = []
    for i in range(min_row, max_row + 1):
        subgrid.append(grid[i][min_col:max_col + 1])

    return subgrid

# Test input grid
input_grid = [
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 5, 2, 0, 0],
    [0, 0, 2, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
]

# Extract and print the subgrid
output_grid = extract_subgrid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))
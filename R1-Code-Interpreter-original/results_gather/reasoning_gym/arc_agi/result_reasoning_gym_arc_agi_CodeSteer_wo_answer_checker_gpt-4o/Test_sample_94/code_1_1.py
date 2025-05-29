def extract_subgrid(input_grid):
    # Determine the dimensions of the input grid
    rows = len(input_grid)
    cols = len(input_grid[0])

    # Initialize variables to find the bounds of the subgrid
    min_row, max_row = rows, 0
    min_col, max_col = cols, 0

    # Iterate over the grid to find the bounds of the non-zero subgrid
    for r in range(rows):
        for c in range(cols):
            if input_grid[r][c] != 0:
                min_row = min(min_row, r)
                max_row = max(max_row, r)
                min_col = min(min_col, c)
                max_col = max(max_col, c)

    # Extract the subgrid using the identified bounds
    subgrid = [input_grid[r][min_col:max_col+1] for r in range(min_row, max_row+1)]

    return subgrid

# Test input grid
input_grid = [
    [0, 6, 5, 6, 6, 9, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 9, 6, 6, 5, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 0, 0, 0, 0],
    [0, 6, 6, 6, 6, 6, 6, 0, 0, 0, 6, 6, 6, 6, 6, 0, 0, 6, 6, 6, 6, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 6, 9, 6, 0, 0, 6, 5, 6, 6, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 5, 6, 6, 6, 0, 0, 6, 6, 6, 6, 0, 0, 0, 0],
    [0, 0, 0, 6, 6, 6, 6, 0, 0, 0, 6, 6, 6, 5, 6, 0, 0, 6, 6, 9, 6, 0, 0, 0, 0],
    [0, 0, 0, 9, 6, 5, 6, 0, 0, 0, 6, 6, 6, 6, 6, 0, 0, 6, 6, 6, 6, 0, 0, 0, 0],
    [0, 0, 0, 6, 6, 9, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 6, 6, 9, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, 6, 6, 6, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 6, 6, 5, 6, 0],
    [6, 6, 6, 6, 6, 6, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 6, 6, 6, 6, 0],
    [6, 6, 9, 6, 6, 6, 0, 0, 0, 6, 6, 6, 5, 6, 9, 6, 6, 0, 0, 0, 6, 9, 6, 6, 0],
    [6, 6, 6, 6, 6, 9, 0, 0, 0, 6, 9, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 6, 5, 6, 6, 6, 0, 0, 0, 6, 6, 6, 5, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0],
    [6, 6, 6, 6, 6, 6, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 5, 6, 6, 9, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0]
]

# Extract and print the subgrid
output_grid = extract_subgrid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))
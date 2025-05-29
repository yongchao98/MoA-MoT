def extract_subgrid(input_grid, rows, cols, subgrid_rows, subgrid_cols):
    # Extract the subgrid from the lower right corner
    subgrid = []
    for i in range(rows - subgrid_rows, rows):
        subgrid.append(input_grid[i][cols - subgrid_cols:cols])
    return subgrid

# Test input grid
input_grid = [
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 2, 2, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 2, 2, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 5, 5, 3, 5, 3, 0, 0, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 5, 5, 3, 5, 3, 0, 0, 3, 5, 3, 2, 2, 3, 5, 3, 5, 5, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 0, 0, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3, 1, 1, 3, 5, 3],
    [3, 0, 0, 3, 5, 3, 5, 5, 3, 5, 3, 5, 5, 3, 5, 3, 1, 1, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3],
    [3, 4, 4, 3, 5, 3, 0, 0, 3, 5, 3, 1, 1, 3, 5, 3, 2, 2, 3, 5, 3],
    [3, 4, 4, 3, 5, 3, 0, 0, 3, 5, 3, 1, 1, 3, 5, 3, 2, 2, 3, 5, 3],
    [3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3, 3, 3, 3, 5, 3]
]

# Dimensions of the input grid
rows = len(input_grid)
cols = len(input_grid[0])

# Based on the examples, let's assume a 3x3 subgrid for this test case
subgrid_rows = 3
subgrid_cols = 3

# Extract the subgrid
output_grid = extract_subgrid(input_grid, rows, cols, subgrid_rows, subgrid_cols)

# Print the output grid
for row in output_grid:
    print(" ".join(map(str, row)))
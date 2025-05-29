def extract_subgrid(grid):
    # Find the column with the distinct pattern (in this case, the column of 2s)
    subgrid = []
    for row in grid:
        if 2 in row:
            # Extract the part of the row starting from the first occurrence of 2
            subgrid.append(row[row.index(2):])
    # Determine the size of the subgrid to match the pattern in the examples
    # In this case, we want the bottom-left corner of the region
    # Let's assume we want the last 4 rows and the first 2 columns of the subgrid
    return [row[:2] for row in subgrid[-4:]]

# Test input grid
input_grid = [
    [1, 9, 1, 1, 1, 1],
    [1, 1, 1, 1, 2, 1],
    [1, 2, 1, 1, 2, 1],
    [1, 1, 1, 1, 2, 1],
    [1, 1, 1, 1, 2, 1],
    [1, 0, 1, 1, 1, 1]
]

# Extract the subgrid
output_grid = extract_subgrid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
def extract_subgrid(input_grid):
    # Initialize bounds for the subgrid
    min_row, max_row = len(input_grid), -1
    min_col, max_col = len(input_grid[0]), -1
    
    # Find the bounds of the subgrid containing all non-9 elements
    for i, row in enumerate(input_grid):
        for j, value in enumerate(row):
            if value != 9:
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
        subgrid.append(input_grid[i][min_col:max_col + 1])
    
    return subgrid

# Test input grid
input_grid = [
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 3, 0, 1, 3, 9],
    [9, 9, 9, 9, 9, 0, 3, 0, 0, 9],
    [9, 0, 0, 3, 9, 1, 0, 1, 3, 9],
    [9, 3, 0, 1, 9, 9, 9, 9, 9, 9],
    [9, 1, 0, 0, 9, 9, 9, 9, 9, 9],
    [9, 0, 3, 1, 9, 9, 0, 0, 0, 0],
    [9, 0, 0, 0, 9, 9, 3, 0, 1, 3],
    [9, 3, 0, 1, 9, 9, 0, 3, 0, 0],
    [9, 9, 9, 9, 9, 9, 0, 0, 0, 1]
]

# Extract and print the output grid
output_grid = extract_subgrid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))
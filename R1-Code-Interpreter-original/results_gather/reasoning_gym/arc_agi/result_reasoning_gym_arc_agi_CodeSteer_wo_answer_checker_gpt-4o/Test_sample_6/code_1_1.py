def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Determine the number of rows and columns
    num_rows = len(input_grid)
    num_cols = len(input_grid[0]) if num_rows > 0 else 0
    
    # Identify the rows and columns that need transformation
    rows_to_transform = set()
    cols_to_transform = set()
    
    for i in range(num_rows):
        for j in range(num_cols):
            if input_grid[i][j] in [6, 9]:
                rows_to_transform.add(i)
                cols_to_transform.add(j)
    
    # Apply the transformation rule to the identified rows and columns
    for i in rows_to_transform:
        for j in range(num_cols):
            if input_grid[i][j] in [6, 9]:
                output_grid[i][j] = 8
    
    for j in cols_to_transform:
        for i in range(num_rows):
            if input_grid[i][j] in [6, 9]:
                output_grid[i][j] = 8
    
    return output_grid

# Test input grid
input_grid = [
    [8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 6, 8],
    [8, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 8],
    [8, 6, 9, 6, 6, 6, 6, 8, 8, 8, 8, 8],
    [8, 6, 6, 6, 6, 9, 6, 8, 6, 6, 6, 6],
    [8, 6, 6, 6, 6, 6, 6, 8, 6, 6, 6, 6],
    [8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 9, 6],
    [8, 8, 8, 8, 6, 6, 6, 8, 6, 6, 6, 6],
    [8, 8, 8, 8, 6, 6, 9, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 6, 6, 6, 8, 8, 8, 8, 8],
    [6, 6, 6, 8, 9, 6, 6, 8, 8, 8, 8, 8],
    [9, 6, 6, 8, 8, 8, 8, 8, 6, 9, 6, 9],
    [6, 6, 6, 8, 8, 8, 8, 8, 6, 6, 6, 6],
    [8, 8, 8, 8, 8, 8, 8, 8, 6, 9, 6, 6]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
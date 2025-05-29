def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Determine the number of rows and columns
    num_rows = len(input_grid)
    num_cols = len(input_grid[0]) if num_rows > 0 else 0
    
    # Function to find contiguous blocks of 6s or 9s
    def find_blocks(grid):
        blocks = []
        for i in range(num_rows):
            start = None
            for j in range(num_cols):
                if grid[i][j] in [6, 9]:
                    if start is None:
                        start = j
                else:
                    if start is not None:
                        blocks.append((i, start, j - 1))
                        start = None
            if start is not None:
                blocks.append((i, start, num_cols - 1))
        return blocks
    
    # Find blocks of 6s or 9s
    blocks = find_blocks(input_grid)
    
    # Apply the transformation rule to the identified blocks
    for i, start_col, end_col in blocks:
        for j in range(start_col, end_col + 1):
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
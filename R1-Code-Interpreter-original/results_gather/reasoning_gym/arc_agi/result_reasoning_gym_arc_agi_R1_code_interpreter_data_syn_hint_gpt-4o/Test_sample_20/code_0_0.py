def transform_grid(input_grid):
    # Initialize the output grid with 2s
    output_grid = [[2 for _ in range(12)] for _ in range(12)]
    
    # Define the blocks based on the input grid
    blocks = [
        (3, 0, 0, 3, 3),  # (value, start_row, start_col, block_height, block_width)
        (9, 1, 1, 3, 3),
        (5, 0, 4, 3, 3)
    ]
    
    # Fill the blocks in the output grid
    for value, start_row, start_col, block_height, block_width in blocks:
        for i in range(start_row, start_row + block_height):
            for j in range(start_col, start_col + block_width):
                output_grid[i][j] = value
    
    return output_grid

# Test input grid
input_grid = [
    [3, 2, 2, 3, 2, 5, 5, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [4, 2, 2, 3, 2, 5, 5, 5]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
def expand_grid(input_grid):
    # Determine the size of the input grid
    input_rows = len(input_grid)
    input_cols = len(input_grid[0])
    
    # Calculate the block size based on the input grid size
    block_size = input_rows
    
    # Calculate the size of the output grid
    output_rows = input_rows * block_size
    output_cols = input_cols * block_size
    
    # Initialize the output grid with 2s
    output_grid = [[2] * output_cols for _ in range(output_rows)]
    
    # Function to place a block in the output grid
    def place_block(value, start_row, start_col):
        for i in range(block_size):
            for j in range(block_size):
                # Ensure we don't go out of bounds
                if 0 <= start_row + i < output_rows and 0 <= start_col + j < output_cols:
                    output_grid[start_row + i][start_col + j] = value
    
    # Iterate over the input grid and place blocks in the output grid
    for i in range(input_rows):
        for j in range(input_cols):
            value = input_grid[i][j]
            if value != 2:
                # Calculate the starting position for the block
                start_row = i * block_size
                start_col = j * block_size
                place_block(value, start_row, start_col)
    
    return output_grid

# Test input grid
input_grid = [
    [3, 2, 2, 3, 2, 5, 5, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [2, 9, 9, 2, 5, 5, 2, 5],
    [4, 2, 2, 3, 2, 5, 5, 5]
]

# Get the output grid
output_grid = expand_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
def transform_grid(input_grid):
    # Find the unique block in the input grid
    unique_block = None
    for i in range(len(input_grid)):
        for j in range(len(input_grid[i])):
            if input_grid[i][j] != 5:  # Assuming 5 is the surrounding number
                if unique_block is None:
                    unique_block = (input_grid[i][j], i, j)
                else:
                    break
        if unique_block is not None:
            break

    # Create the output grid by replicating the unique block in the middle
    output_grid = [row[:] for row in input_grid]
    block_value, block_row, block_col = unique_block
    for i in range(len(input_grid)):
        for j in range(len(input_grid[i])):
            if input_grid[i][j] == block_value:
                output_grid[i][j] = block_value

    # Replicate the block in the middle
    mid_row = len(input_grid) // 2
    mid_col = len(input_grid[0]) // 2
    for i in range(2):
        for j in range(2):
            output_grid[mid_row - 1 + i][mid_col - 1 + j] = block_value

    return output_grid

# Test input grid
input_grid = [
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 0, 0, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
]

output_grid = transform_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))
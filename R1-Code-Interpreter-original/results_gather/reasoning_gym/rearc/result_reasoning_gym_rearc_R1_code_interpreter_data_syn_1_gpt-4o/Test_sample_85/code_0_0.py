def transform_grid(input_grid):
    # Copy the input grid to the output grid
    output_grid = [row[:] for row in input_grid]
    
    # Define the blocks and their middle positions
    blocks = [
        ((3, 7), (3, 8), (3, 9)),
        ((4, 7), (4, 8), (4, 9)),
        ((5, 7), (5, 8), (5, 9))
    ]
    
    # Change the middle element of each block
    for block in blocks:
        middle = block[1]  # Middle element of the block
        output_grid[middle[0]][middle[1]] = 1  # Change to the surrounding number
    
    return output_grid

# Test input grid
input_grid = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 1],
    [1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 1],
    [1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
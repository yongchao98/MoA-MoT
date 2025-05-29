def transform_grid(input_grid):
    # Initialize the output grid
    output_grid = [[8] * len(input_grid[0]) for _ in range(3)]
    
    # Find the position for '1' and '3' in the first and third rows
    for row in [0, 2]:
        found_nine = False
        for i in range(len(input_grid[0])):
            if input_grid[0][i] == 9 and not found_nine:
                output_grid[row][i] = 1
                found_nine = True
            elif found_nine and input_grid[0][i] == 8:
                output_grid[row][i] = 3
                break
    
    return output_grid

# Test input grid
input_grid = [
    [9, 9, 9, 8, 8, 8, 8, 8, 8],
    [8, 9, 8, 8, 8, 8, 8, 9, 8],
    [8, 8, 8, 9, 8, 8, 9, 9, 9],
    [8, 8, 9, 9, 9, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 9, 8, 8, 8],
    [9, 9, 9, 8, 9, 9, 9, 8, 8],
    [8, 9, 8, 8, 8, 8, 8, 8, 8],
    [8, 8, 8, 8, 8, 8, 9, 9, 9],
    [8, 8, 9, 9, 9, 8, 8, 9, 8],
    [8, 8, 8, 9, 8, 8, 8, 8, 8]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
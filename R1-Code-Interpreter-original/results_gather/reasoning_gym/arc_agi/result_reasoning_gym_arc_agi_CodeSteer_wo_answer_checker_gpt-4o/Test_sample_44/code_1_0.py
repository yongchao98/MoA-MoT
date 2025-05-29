def transform_grid(input_grid):
    # Find the largest number in the grid
    largest_number = max(max(row) for row in input_grid)
    
    # Initialize the output grid with the same dimensions
    output_grid = [[0] * len(input_grid[0]) for _ in range(len(input_grid))]
    
    # Place the largest number in the first row and last column
    for i in range(len(input_grid)):
        for j in range(len(input_grid[0])):
            if i == 0 or j == len(input_grid[0]) - 1:
                output_grid[i][j] = largest_number
            else:
                # Fill the remaining positions with the corresponding input values
                output_grid[i][j] = input_grid[i][j]
    
    return output_grid

# Test input grid
input_grid = [
    [1, 8, 8],
    [1, 8, 8],
    [8, 1, 8]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
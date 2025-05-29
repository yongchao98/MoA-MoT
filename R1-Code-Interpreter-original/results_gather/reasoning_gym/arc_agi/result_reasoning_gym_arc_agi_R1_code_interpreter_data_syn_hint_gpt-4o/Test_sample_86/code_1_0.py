def transform_grid(input_grid):
    # Initialize the output grid
    output_grid = [[8] * len(input_grid[0]) for _ in range(3)]
    
    # Function to find the positions for '1' and '3'
    def find_positions(row):
        found_nine = False
        for i in range(len(row)):
            if row[i] == 9 and not found_nine:
                output_grid[0][i] = 1
                output_grid[2][i] = 1
                found_nine = True
            elif found_nine and row[i] == 8:
                output_grid[0][i] = 3
                output_grid[2][i] = 3
                break
    
    # Apply the function to the first row of the input grid
    find_positions(input_grid[0])
    
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
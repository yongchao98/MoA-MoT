def transform_grid(input_grid):
    # Determine the number of rows and columns
    num_rows = len(input_grid)
    num_cols = len(input_grid[0])
    
    # Initialize the output grid
    output_grid = [[0] * num_cols for _ in range(num_rows)]
    
    # Transform the last column to a sequence starting from 1
    for row in range(num_rows):
        output_grid[row][-1] = 1
    
    # Transform the third column to 2 if it contains specific values
    for row in range(num_rows):
        if input_grid[row][2] in [5, 7]:
            output_grid[row][2] = 2
        else:
            output_grid[row][2] = input_grid[row][2]
    
    # Transform the first column to 3 if it contains specific values
    for row in range(num_rows):
        if input_grid[row][0] in [1, 5]:
            output_grid[row][0] = 3
        else:
            output_grid[row][0] = input_grid[row][0]
    
    # Copy the unchanged columns
    for row in range(num_rows):
        for col in range(num_cols):
            if col not in [0, 2, num_cols - 1]:
                output_grid[row][col] = input_grid[row][col]
    
    return output_grid

# Test input grid
input_grid = [
    [5, 8, 5, 5, 5],
    [5, 8, 5, 5, 5],
    [5, 8, 5, 5, 9],
    [5, 8, 7, 5, 9],
    [1, 8, 7, 5, 9]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
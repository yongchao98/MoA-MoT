def transform_grid(input_grid):
    # Define the dimensions of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Initialize the output grid with the same dimensions
    output_grid = [[0] * cols for _ in range(rows)]
    
    # Apply the transformation rules
    for i in range(rows):
        for j in range(cols):
            # Copy the input value to the output grid
            output_grid[i][j] = input_grid[i][j]
            
            # Apply specific transformations based on the input value
            if input_grid[i][j] == 4:
                # Example transformation for 4
                if j + 1 < cols:
                    output_grid[i][j + 1] = 9
            elif input_grid[i][j] == 5:
                # Example transformation for 5
                if i + 1 < rows:
                    output_grid[i + 1][j] = 8
            elif input_grid[i][j] == 9:
                # Example transformation for 9
                if i - 1 >= 0:
                    output_grid[i - 1][j] = 2
    
    return output_grid

# Test input grid
input_grid = [
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 9, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [2, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 9, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4],
    [4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 5, 5, 5, 5, 4, 4]
]

# Transform the input grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
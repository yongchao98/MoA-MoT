def transform_grid(input_grid):
    # Get the dimensions of the grid
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Create a copy of the input grid to modify
    output_grid = [row[:] for row in input_grid]
    
    # Iterate over the grid to find special numbers
    for i in range(rows):
        for j in range(cols):
            # Check if the current number is a special number
            if input_grid[i][j] != input_grid[0][0]:
                special_number = input_grid[i][j]
                
                # Apply the transformation pattern
                # Center
                output_grid[i][j] = special_number
                
                # Up
                if i > 0:
                    output_grid[i-1][j] = 1
                
                # Down
                if i < rows - 1:
                    output_grid[i+1][j] = 1
                
                # Left
                if j > 0:
                    output_grid[i][j-1] = 1
                
                # Right
                if j < cols - 1:
                    output_grid[i][j+1] = 1
                
                # Top-left
                if i > 0 and j > 0:
                    output_grid[i-1][j-1] = special_number
                
                # Top-right
                if i > 0 and j < cols - 1:
                    output_grid[i-1][j+1] = special_number
                
                # Bottom-left
                if i < rows - 1 and j > 0:
                    output_grid[i+1][j-1] = special_number
                
                # Bottom-right
                if i < rows - 1 and j < cols - 1:
                    output_grid[i+1][j+1] = special_number
    
    return output_grid

# Test input grid
input_grid = [
    [3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 8, 3],
    [3, 3, 3, 8, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3],
    [3, 3, 3, 3, 3, 3, 3, 3]
]

# Transform the grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
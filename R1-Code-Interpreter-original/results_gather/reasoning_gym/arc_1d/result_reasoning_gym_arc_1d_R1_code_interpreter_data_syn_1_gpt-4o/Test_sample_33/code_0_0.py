def transform_grid(input_grid):
    # Assuming the grid is 5x4 based on the input length
    rows = 5
    cols = 4
    output_grid = input_grid[:]
    
    # Fill the second and third rows with the first row if the first row has non-zero numbers
    if any(input_grid[i] != 0 for i in range(cols)):
        for row in range(1, 3):
            for col in range(cols):
                output_grid[row * cols + col] = input_grid[col]
    
    return output_grid

# Test input
input_grid = [2, 2, 2, 2, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0]
output_grid = transform_grid(input_grid)
print(output_grid)
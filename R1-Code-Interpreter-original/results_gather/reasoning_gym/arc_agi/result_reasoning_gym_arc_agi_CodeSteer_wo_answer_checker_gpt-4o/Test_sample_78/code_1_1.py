def generate_output_grid(input_grid):
    # Determine the size of the output grid
    size = len(input_grid)
    
    # Create the output grid filled with 9s
    output_grid = [[9] * size for _ in range(size)]
    
    # Fill the diagonal with the input numbers
    for i in range(size):
        output_grid[i][i] = input_grid[i]
    
    return output_grid

# Test input
input_grid = [6, 7, 4, 5, 9]
output_grid = generate_output_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
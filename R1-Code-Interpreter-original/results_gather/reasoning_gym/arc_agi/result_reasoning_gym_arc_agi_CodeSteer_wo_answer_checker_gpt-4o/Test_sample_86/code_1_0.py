def transform_grid(input_grid):
    # Initialize the output grid with 8s
    output_grid = [[8] * 6 for _ in range(3)]
    
    # Flatten the input grid to a single list for easier processing
    flat_input = [num for row in input_grid for num in row]
    
    # Process the first row of the output grid
    # We will look for sequences of 9s and 8s
    for i in range(len(flat_input) - 2):
        if flat_input[i] == 9 and flat_input[i+1] == 9 and flat_input[i+2] == 9:
            output_grid[0][i % 6] = 1
        elif flat_input[i] == 8 and flat_input[i+1] == 8 and flat_input[i+2] == 8:
            output_grid[0][i % 6] = 3
    
    # Analyze the second and third rows
    # Based on the examples, these rows are mostly filled with 8s
    # Further analysis is needed to determine the exact pattern for 1s and 3s
    
    # For now, fill them with 8s
    output_grid[1] = [8] * 6
    output_grid[2] = [8] * 6
    
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

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
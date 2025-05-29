def extract_output_grid(input_grid):
    # Define the size of the output grid based on the examples
    output_height = 11
    output_width = 8
    
    # Initialize the output grid with 7s
    output_grid = [[7] * output_width for _ in range(output_height)]
    
    # Analyze the input grid to fill the output grid
    # This is a hypothesis based on the examples
    for i in range(output_height):
        for j in range(output_width):
            # Example pattern: take the minimum of a 2x2 block from the input grid
            # This is a guess based on the examples
            block = [input_grid[i*2][j*2], input_grid[i*2][j*2+1],
                     input_grid[i*2+1][j*2], input_grid[i*2+1][j*2+1]]
            output_grid[i][j] = min(block)
    
    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 6, 7, 7, 7, 7, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 7],
    [7, 7, 1, 7, 7, 7, 7, 4, 2, 7, 4, 4, 7, 4, 9, 7, 4, 9, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 7],
    [7, 2, 7, 7, 7, 7, 7, 4, 2, 7, 4, 6, 7, 4, 6, 7, 4, 4, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 1, 7, 7, 7, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 4, 2, 7, 4, 6, 7, 4, 4, 7, 4, 4, 7, 7],
    [7, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 4, 4, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 4, 0, 7, 4, 4, 7, 4, 4, 7, 4, 0, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 9, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 1, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
]

output_grid = extract_output_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
def transform_grid(input_grid):
    # Initialize the output grid with the same dimensions as the input grid
    output_grid = [[7 for _ in range(len(input_grid[0]))] for _ in range(len(input_grid))]
    
    # Step 1: Expand 9s into 2x2 blocks
    for i in range(len(input_grid)):
        for j in range(len(input_grid[0])):
            if input_grid[i][j] == 9:
                for di in range(2):
                    for dj in range(2):
                        if 0 <= i + di < len(input_grid) and 0 <= j + dj < len(input_grid[0]):
                            output_grid[i + di][j + dj] = 9
    
    # Step 2: Introduce 3s in a pattern
    for i in range(len(input_grid)):
        for j in range(len(input_grid[0])):
            if output_grid[i][j] == 7:  # Only place 3s where there are no 9s
                if (i % 4 == 0 or j % 4 == 0):
                    output_grid[i][j] = 3
    
    return output_grid

# Test input grid
test_input = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 9, 9, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
]

# Transform the test input grid
output_grid = transform_grid(test_input)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
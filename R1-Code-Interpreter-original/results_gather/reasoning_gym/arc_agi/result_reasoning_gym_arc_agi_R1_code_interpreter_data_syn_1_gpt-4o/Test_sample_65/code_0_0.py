def generate_output_grid(input_grid):
    # Find the primary number in the last row
    primary_number = None
    for num in input_grid[-1]:
        if num != 1:
            primary_number = num
            break
    
    # Initialize the output grid
    output_grid = [[1] * 10 for _ in range(10)]
    
    # Fill the output grid based on the rule
    for i in range(10):
        for j in range(10):
            if (i + j) % 2 == 0:
                output_grid[i][j] = primary_number
    
    # Place '5' in the specific positions
    output_grid[0][0] = 5
    output_grid[0][9] = 5
    output_grid[9][0] = 5
    output_grid[9][9] = 5
    
    return output_grid

# Test input grid
input_grid = [
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1, 1, 1, 3, 1, 1]
]

output_grid = generate_output_grid(input_grid)
for row in output_grid:
    print(' '.join(map(str, row)))
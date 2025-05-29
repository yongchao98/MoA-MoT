def transform_grid(input_grid):
    # Define the positions to check in the input grid for the 3x3 output grid
    positions = [
        (0, 0), (0, 5), (0, 10),
        (5, 0), (5, 5), (5, 10),
        (10, 0), (10, 5), (10, 10)
    ]
    
    # Initialize the output grid
    output_grid = []
    
    # Check each position and determine the output value
    for i in range(3):
        row = []
        for j in range(3):
            x, y = positions[i * 3 + j]
            if input_grid[x][y] == 5:
                row.append(6)
            else:
                row.append(3)
        output_grid.append(row)
    
    return output_grid

# Test input grid
input_grid = [
    [3, 3, 3, 2, 3, 5, 3, 2, 3, 3, 3],
    [5, 3, 3, 2, 3, 3, 3, 2, 3, 5, 3],
    [3, 3, 5, 2, 3, 3, 3, 2, 3, 3, 3],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [3, 3, 3, 2, 5, 3, 3, 2, 3, 3, 5],
    [3, 3, 5, 2, 3, 5, 3, 2, 3, 3, 3],
    [3, 3, 3, 2, 3, 3, 3, 2, 3, 5, 3],
    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3],
    [3, 5, 3, 2, 3, 3, 5, 2, 3, 3, 3],
    [3, 3, 3, 2, 3, 3, 3, 2, 5, 3, 5]
]

# Get the output grid
output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
# Define the input grid
input_grid = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0, 0, 2, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0]
]

# Initialize the output grid with zeros
output_grid = [[0 for _ in range(9)] for _ in range(9)]

# Apply the rule: retain non-zero blocks or distinct numbers
for i in range(9):
    for j in range(9):
        if input_grid[i][j] != 0:
            output_grid[i][j] = input_grid[i][j]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
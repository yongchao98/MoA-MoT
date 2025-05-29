# Test input grid
input_grid = [
    [2, 0, 7, 7, 0, 1],
    [0, 0, 0, 0, 0, 0],
    [7, 0, 3, 3, 0, 7],
    [7, 0, 7, 7, 0, 7],
    [7, 0, 3, 3, 0, 7],
    [7, 0, 7, 7, 0, 7],
    [0, 0, 0, 0, 0, 0],
    [9, 0, 7, 7, 0, 8]
]

# Function to generate the output grid based on the identified pattern
def generate_output_grid(input_grid):
    output_grid = []
    for i in range(0, len(input_grid), 2):
        row = []
        for j in range(0, len(input_grid[0]), 3):
            row.append(input_grid[i][j])
        output_grid.append(row)
    return output_grid

# Generate and print the output grid
output_grid = generate_output_grid(input_grid)
for row in output_grid:
    print(" ".join(map(str, row)))
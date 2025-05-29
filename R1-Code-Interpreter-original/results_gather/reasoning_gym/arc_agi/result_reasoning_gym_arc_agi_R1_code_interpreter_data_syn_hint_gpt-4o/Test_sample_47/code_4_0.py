# Define the input grid
input_grid = [
    [6, 6, 6, 9, 6, 6, 6, 9, 6, 1, 6],
    [6, 1, 6, 9, 6, 6, 6, 9, 6, 6, 6],
    [6, 6, 6, 9, 6, 1, 6, 9, 1, 6, 6],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [6, 6, 6, 9, 6, 6, 1, 9, 6, 6, 6],
    [6, 6, 6, 9, 6, 1, 6, 9, 6, 6, 1],
    [6, 1, 6, 9, 6, 6, 6, 9, 6, 6, 6],
    [9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9],
    [1, 6, 6, 9, 6, 6, 6, 9, 6, 6, 6],
    [6, 6, 6, 9, 1, 6, 6, 9, 6, 1, 6],
    [1, 6, 6, 9, 6, 6, 1, 9, 6, 6, 6]
]

# Initialize the output grid
output_grid = [[6, 6, 6], [6, 6, 6], [6, 6, 6]]

# Define the regions to check for '1' in the input grid
regions = [
    [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)],  # Top-left 3x3 section
    [(4, 4), (4, 5), (4, 6), (5, 4), (5, 5), (5, 6), (6, 4), (6, 5), (6, 6)],  # Middle 3x3 section
    [(8, 8), (8, 9), (8, 10), (9, 8), (9, 9), (9, 10), (10, 8), (10, 9), (10, 10)] # Bottom-right 3x3 section
]

# Update the output grid based on the presence of '1' in the input grid
for i, region in enumerate(regions):
    for pos in region:
        if input_grid[pos[0]][pos[1]] == 1:
            output_grid[i // 3][i % 3] = 4

# Print the output grid
print(output_grid)
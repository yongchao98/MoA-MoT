# Define the pattern of 7s identified in the input grid
pattern = [
    [7, 7, 7],
    [7, 9, 7],
    [9, 7, 9]
]

# Repeat the pattern to form the output grid
output_grid = []
for _ in range(3):  # Repeat vertically
    for row in pattern:
        output_grid.append(row)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
# Define the input grid
input_grid = [
    [4],
    [4]
]

# Initialize the output grid
output_grid = []

# Apply the rule to generate the output grid
for row in input_grid:
    # Repeat each row twice
    new_row = row * 2
    output_grid.append(new_row)
    output_grid.append(new_row)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
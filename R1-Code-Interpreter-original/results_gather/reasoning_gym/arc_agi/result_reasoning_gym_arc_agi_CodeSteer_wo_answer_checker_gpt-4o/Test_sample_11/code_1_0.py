# Define the input grid
input_grid = [
    [4, 8, 8, 8, 8, 8, 8, 8, 8, 6],
    [4, 8, 0, 8, 8, 8, 0, 8, 8, 6],
    [4, 8, 8, 8, 0, 8, 8, 8, 8, 6],
    [4, 8, 8, 8, 8, 8, 8, 8, 8, 6],
    [4, 8, 8, 0, 8, 8, 8, 8, 8, 6],
    [4, 8, 8, 8, 8, 8, 0, 8, 8, 6],
    [4, 8, 8, 0, 8, 8, 8, 8, 8, 6],
    [4, 8, 8, 8, 8, 8, 8, 8, 8, 6],
    [4, 8, 8, 8, 8, 0, 8, 8, 0, 6],
    [4, 0, 8, 8, 8, 8, 8, 8, 8, 6]
]

# Get the replacement number from the first element of the first row
replacement_number = input_grid[0][0]

# Create the output grid by applying the rule
output_grid = []
for i, row in enumerate(input_grid):
    if i == len(input_grid) - 1:
        # Keep the last row unchanged
        output_grid.append(row)
    else:
        # Replace '0' with the replacement number
        new_row = [replacement_number if x == 0 else x for x in row]
        output_grid.append(new_row)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
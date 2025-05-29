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

# Apply the rule to generate the output grid
output_grid = []
for row in input_grid:
    first_number = row[0]
    new_row = [first_number if x == 0 else x for x in row]
    output_grid.append(new_row)

# Print the output grid
print(output_grid)
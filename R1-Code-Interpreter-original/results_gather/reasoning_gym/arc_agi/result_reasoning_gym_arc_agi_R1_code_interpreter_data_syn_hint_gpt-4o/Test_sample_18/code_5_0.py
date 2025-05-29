def construct_output_grid(input_grid):
    # Initialize an empty output grid
    output_grid = []

    # Iterate over each row in the input grid
    for row in input_grid:
        # Create a new row for the output grid
        new_row = []
        for value in row:
            # If the value is 7, add it to the new row
            if value == 7:
                new_row.append(7)
            else:
                # Otherwise, add 9
                new_row.append(9)
        # Add the new row to the output grid
        output_grid.append(new_row)

    return output_grid

# Define the test input grid
input_grid = [
    [9, 9, 9, 9, 2, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 2, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9],
    [7, 7, 7, 9, 9, 9, 9, 9, 9],
    [7, 9, 7, 9, 9, 9, 9, 9, 9],
    [9, 7, 9, 9, 9, 9, 9, 9, 9],
    [9, 9, 9, 9, 9, 9, 9, 9, 9]
]

# Construct the output grid
output_grid = construct_output_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
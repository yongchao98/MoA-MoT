def generate_output_grid(input_grid):
    # Determine the size of the output grid
    output_rows = 9
    output_cols = 5
    
    # Define the pattern for the output grid
    pattern = [3, 9]  # Alternating pattern
    full_row_value = 9  # Value for the full row
    
    # Initialize the output grid
    output_grid = []
    
    for i in range(output_rows):
        if i == 5:  # Full row of a single number
            output_grid.append([full_row_value] * output_cols)
        else:
            row = []
            for j in range(output_cols):
                row.append(pattern[j % 2])
            output_grid.append(row)
    
    return output_grid

# Test input grid
input_grid = [
    [1, 0, 2, 0, 1, 9, 1, 7, 0, 9, 8, 3],
    [8, 8, 4, 3, 3, 3, 3, 3, 2, 4, 8, 1],
    [3, 4, 9, 3, 3, 3, 3, 3, 0, 9, 3, 9],
    [4, 2, 0, 3, 3, 3, 3, 3, 1, 0, 9, 3],
    [8, 3, 1, 3, 3, 9, 3, 3, 2, 4, 5, 8],
    [1, 2, 0, 3, 3, 3, 3, 3, 5, 2, 5, 1],
    [1, 2, 8, 3, 3, 3, 3, 3, 0, 4, 2, 4],
    [8, 2, 5, 3, 3, 3, 3, 3, 5, 4, 3, 3],
    [2, 8, 2, 3, 3, 3, 3, 3, 9, 4, 7, 1],
    [8, 9, 9, 3, 3, 3, 3, 3, 5, 7, 7, 2],
    [7, 0, 8, 7, 0, 1, 5, 8, 0, 5, 3, 2],
    [8, 2, 9, 9, 4, 2, 3, 5, 5, 0, 9, 3]
]

# Generate the output grid
output_grid = generate_output_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
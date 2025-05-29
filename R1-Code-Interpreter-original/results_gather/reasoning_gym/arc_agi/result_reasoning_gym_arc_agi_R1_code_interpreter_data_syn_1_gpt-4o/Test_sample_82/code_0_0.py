# Define the unique numbers forming clusters in the input grid
unique_numbers = [0, 3, 6, 8, 9]

# Create the output grid by repeating the sequence of unique numbers
output_grid = [unique_numbers for _ in range(10)]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
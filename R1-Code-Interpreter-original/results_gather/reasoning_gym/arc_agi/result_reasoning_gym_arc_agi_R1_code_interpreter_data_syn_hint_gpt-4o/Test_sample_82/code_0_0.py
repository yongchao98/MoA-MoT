# Define the unique numbers found in the test input grid
unique_numbers = [3, 6, 8, 9]

# Sort the unique numbers
unique_numbers.sort()

# Create the output grid by repeating the sorted unique numbers
output_grid = [unique_numbers for _ in range(len(unique_numbers))]

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))
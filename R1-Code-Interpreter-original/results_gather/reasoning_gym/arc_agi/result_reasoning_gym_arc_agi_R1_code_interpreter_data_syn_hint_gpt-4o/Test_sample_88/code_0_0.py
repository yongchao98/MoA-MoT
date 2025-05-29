# Define the unique numbers from the test input grid
unique_numbers = [3, 6]

# Generate the output grid by forming pairs of the unique numbers
output_grid = [(unique_numbers[i], unique_numbers[j]) for i in range(len(unique_numbers)) for j in range(i, len(unique_numbers))]

# Print the output grid
print(output_grid)
# Define the input grid
input_grid = [0, 7, 7, 7, 7, 0, 0, 0, 0, 7, 7, 7, 7, 7, 0, 0, 0, 0, 7, 7, 0, 0, 0, 0]

# Extract non-zero elements
non_zero_elements = [x for x in input_grid if x != 0]

# Calculate the number of zeros needed
num_zeros = len(input_grid) - len(non_zero_elements)

# Create the output grid
output_grid = non_zero_elements + [0] * num_zeros

# Print the output grid
print(output_grid)
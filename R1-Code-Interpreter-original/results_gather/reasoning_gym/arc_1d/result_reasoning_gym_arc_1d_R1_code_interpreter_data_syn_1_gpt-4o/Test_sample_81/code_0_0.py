# Define the test input
test_input = [0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Find the first non-zero element and its position
start_position = next((i for i, x in enumerate(test_input) if x != 0), None)
sequence_length = sum(1 for x in test_input if x != 0)
grid_length = len(test_input)

# Calculate the shift
shift = grid_length - (start_position + sequence_length)

# Create the output grid
output_grid = [0] * grid_length
for i in range(sequence_length):
    output_grid[start_position + shift + i] = test_input[start_position + i]

# Print the output grid
print(output_grid)
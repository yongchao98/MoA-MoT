# Define the original 5x5 matrix
matrix = [
    [-9, 20, -4, 6, -9],
    [-2, 10, 14, 19, 0],
    [-3, 12, 15, 7, -10],
    [6, -9, 17, 7, 11],
    [8, 2, 6, -9, -5]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = [[0] * output_size for _ in range(output_size)]

# Perform average pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the 3x3 window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Flatten the window and calculate the average
        window_values = [value for row in window for value in row]
        average = sum(window_values) / len(window_values)
        # Store the average in the output matrix
        output_matrix[i][j] = average

# Print the output matrix
print(output_matrix)
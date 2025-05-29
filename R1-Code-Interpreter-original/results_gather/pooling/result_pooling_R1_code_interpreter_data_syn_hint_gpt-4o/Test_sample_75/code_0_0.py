# Define the original 5x5 matrix
matrix = [
    [6, 6, -10, -3, -7],
    [-7, 15, -7, 4, 20],
    [12, -3, -4, 5, 12],
    [17, 12, 20, -6, -5],
    [3, 0, -3, -6, 17]
]

# Define the size of the sliding window
window_size = 3

# Calculate the dimensions of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = [[0] * output_size for _ in range(output_size)]

# Perform the average pooling operation
for i in range(output_size):
    for j in range(output_size):
        # Extract the 3x3 window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Flatten the window and calculate the average
        window_values = [num for row in window for num in row]
        average = sum(window_values) / len(window_values)
        # Store the average in the output matrix
        output_matrix[i][j] = average

# Print the output matrix
print(output_matrix)
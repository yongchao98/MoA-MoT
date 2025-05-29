# Define the original matrix
matrix = [
    [10, -5, 20, 11, 6],
    [5, 18, 7, 1, 16],
    [17, -2, 1, -2, 11],
    [7, 17, 14, 17, 9],
    [17, 15, 4, -1, -8]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = [[0] * output_size for _ in range(output_size)]

# Perform the average pooling operation
for i in range(output_size):
    for j in range(output_size):
        # Extract the current window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Flatten the window and calculate the average
        window_values = [num for row in window for num in row]
        average = sum(window_values) / len(window_values)
        # Store the average in the output matrix
        output_matrix[i][j] = average

# Print the output matrix
print(output_matrix)
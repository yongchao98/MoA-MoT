# Define the original matrix
matrix = [
    [0, 12, 3, 11, -7],
    [7, 10, -2, -7, -3],
    [15, 18, 5, 14, 20],
    [-1, 6, 4, 16, 8],
    [-8, 14, -5, 10, 3]
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
        # Extract the current window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Flatten the window and calculate the average
        window_values = [num for row in window for num in row]
        average = sum(window_values) / len(window_values)
        # Store the average in the output matrix
        output_matrix[i][j] = average

# Print the output matrix
print(output_matrix)
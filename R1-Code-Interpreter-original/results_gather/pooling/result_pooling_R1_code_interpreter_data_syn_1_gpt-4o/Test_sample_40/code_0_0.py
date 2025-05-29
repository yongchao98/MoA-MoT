# Define the original matrix
matrix = [
    [-4, -9, -1, 7, -8],
    [13, 0, 8, 8, 9],
    [0, 19, -6, 0, 19],
    [7, 12, 12, -1, -3],
    [9, -5, 5, 11, 0]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the resulting matrix
result_size = len(matrix) - window_size + 1

# Initialize the result matrix
result_matrix = [[0] * result_size for _ in range(result_size)]

# Perform max pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the current window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Find the maximum value in the current window
        max_value = max(max(row) for row in window)
        # Store the maximum value in the result matrix
        result_matrix[i][j] = max_value

# Print the result matrix
print(result_matrix)
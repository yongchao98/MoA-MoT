# Define the original 5x5 matrix
matrix = [
    [2, 12, 10, 5, 7],
    [19, -5, 17, 20, -1],
    [2, 20, 18, 20, 2],
    [0, -9, -1, -3, -5],
    [10, 5, -6, -5, -2]
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
        # Extract the current 3x3 window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Find the maximum value in the current window
        max_value = max(max(row) for row in window)
        # Store the maximum value in the result matrix
        result_matrix[i][j] = max_value

# Print the result matrix
print(result_matrix)
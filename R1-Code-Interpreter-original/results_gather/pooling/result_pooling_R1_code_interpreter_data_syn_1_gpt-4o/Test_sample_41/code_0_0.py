# Define the original matrix
matrix = [
    [6, 3, 8, 4, 20],
    [-8, -8, 0, 19, -10],
    [-7, 15, 2, -2, 10],
    [9, -9, 10, 7, 13],
    [-5, 11, -9, -1, 12]
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
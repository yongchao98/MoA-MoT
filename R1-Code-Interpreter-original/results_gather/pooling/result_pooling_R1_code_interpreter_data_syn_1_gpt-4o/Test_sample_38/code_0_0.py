# Define the original matrix
matrix = [
    [0, 9, -1, -10, 20],
    [4, 2, -7, -3, 6],
    [-5, -10, 5, 14, 3],
    [7, -4, 7, -1, 6],
    [7, 5, 2, -9, 17]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = [[0] * output_size for _ in range(output_size)]

# Perform max pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the current window
        current_window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Find the maximum value in the current window
        max_value = max(max(row) for row in current_window)
        # Store the maximum value in the output matrix
        output_matrix[i][j] = max_value

# Print the output matrix
print(output_matrix)
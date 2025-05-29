# Define the original matrix
matrix = [
    [-7, 19, -9, -8, 6],
    [6, 11, 4, -6, -10],
    [0, 1, 6, 16, -9],
    [-10, 13, 14, -8, -1],
    [12, 18, 10, 20, 5]
]

# Define the size of the sliding window
window_size = 3
n = len(matrix)  # Size of the original matrix

# Initialize the result matrix
result_matrix = []

# Perform max pooling
for i in range(n - window_size + 1):
    row = []
    for j in range(n - window_size + 1):
        # Extract the current 3x3 window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Find the maximum value in the current window
        max_value = max(max(row) for row in window)
        row.append(max_value)
    result_matrix.append(row)

# Print the result matrix
print(result_matrix)
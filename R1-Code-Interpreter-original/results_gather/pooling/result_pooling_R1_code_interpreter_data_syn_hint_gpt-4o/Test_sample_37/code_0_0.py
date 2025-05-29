# Define the original matrix
matrix = [
    [0, 19, -7, 4, -4],
    [-1, 14, 2, 5, 18],
    [10, -7, 11, 13, 16],
    [20, 19, 1, 6, -7],
    [11, -9, -3, 1, -3]
]

# Define the size of the sliding window
window_size = 3
n = len(matrix)

# Initialize the result matrix
result = []

# Perform max pooling
for i in range(n - window_size + 1):
    row = []
    for j in range(n - window_size + 1):
        # Extract the 3x3 window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Find the maximum value in the window
        max_value = max(max(row) for row in window)
        row.append(max_value)
    result.append(row)

# Print the result matrix
print(result)
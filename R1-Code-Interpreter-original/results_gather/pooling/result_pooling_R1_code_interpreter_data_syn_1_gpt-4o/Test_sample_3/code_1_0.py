# Define the matrix
matrix = [
    [1, 5, 8],
    [6, 15, -2],
    [-3, -4, -3]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the resulting matrix
result_size = len(matrix) - window_size + 1

# Initialize the resulting matrix
result = [[0] * result_size for _ in range(result_size)]

# Perform max pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the current window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Find the maximum value in the current window
        result[i][j] = max(window)

# Print the resulting matrix
print(result)
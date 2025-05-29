# Define the matrix
matrix = [
    [-6, 12, 13],
    [-5, -3, 8],
    [4, 8, 18]
]

# Define the size of the sliding window
window_size = 2

# Initialize the result matrix
result = []

# Perform max pooling
for i in range(len(matrix) - window_size + 1):
    row = []
    for j in range(len(matrix[0]) - window_size + 1):
        # Extract the current window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Find the maximum value in the current window
        max_value = max(window)
        row.append(max_value)
    result.append(row)

# Print the result matrix
print(result)
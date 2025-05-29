# Define the original matrix
matrix = [
    [1, 3, 8],
    [13, -7, 4],
    [-1, -9, -4]
]

# Define the size of the sliding window
window_size = 2

# Initialize the resulting matrix
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

# Print the resulting matrix
print(result)
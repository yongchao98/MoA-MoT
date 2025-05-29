# Define the matrix
matrix = [
    [19, -6, 14],
    [-10, 8, -10],
    [12, 1, -10]
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
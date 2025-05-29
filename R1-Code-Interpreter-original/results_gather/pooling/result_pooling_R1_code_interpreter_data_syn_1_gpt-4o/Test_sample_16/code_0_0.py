# Define the matrix
matrix = [
    [10, 18, 10],
    [19, -6, -10],
    [1, -9, 18]
]

# Define the size of the sliding window
window_size = 2

# Initialize the result matrix
result = []

# Perform max pooling
for i in range(len(matrix) - window_size + 1):
    row = []
    for j in range(len(matrix[0]) - window_size + 1):
        # Extract the 2x2 window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Find the maximum value in the window
        max_value = max(window)
        row.append(max_value)
    result.append(row)

# Print the result matrix
print(result)
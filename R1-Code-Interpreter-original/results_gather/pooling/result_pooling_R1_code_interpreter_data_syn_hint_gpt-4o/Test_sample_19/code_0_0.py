# Define the original matrix
matrix = [
    [4, -2, -2, -3],
    [-8, 17, 1, 19],
    [15, 14, -9, 7],
    [-9, 3, 3, 19]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the new matrix
new_size = len(matrix) - window_size + 1

# Initialize the new matrix
new_matrix = [[0] * new_size for _ in range(new_size)]

# Perform max pooling
for i in range(new_size):
    for j in range(new_size):
        # Extract the current 2x2 window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Find the maximum value in the current window
        new_matrix[i][j] = max(window)

# Print the resulting matrix
print(new_matrix)
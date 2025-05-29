# Define the original matrix
matrix = [
    [5, 3, 2, 7],
    [-9, 13, 10, -3],
    [12, 3, -2, 16],
    [14, 4, 9, 18]
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

# Print the new matrix
print(new_matrix)
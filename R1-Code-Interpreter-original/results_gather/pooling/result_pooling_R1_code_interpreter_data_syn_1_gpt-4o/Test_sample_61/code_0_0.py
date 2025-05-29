# Define the original matrix
matrix = [
    [-1, 11, -10, 8],
    [8, 13, -9, 13],
    [-8, 12, 5, 0],
    [-6, -6, 20, -9]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the new matrix
new_size = len(matrix) - window_size + 1

# Initialize the new matrix to store the averages
new_matrix = [[0] * new_size for _ in range(new_size)]

# Perform average pooling
for i in range(new_size):
    for j in range(new_size):
        # Extract the 2x2 window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Calculate the average of the window
        average = sum(window) / len(window)
        # Store the average in the new matrix
        new_matrix[i][j] = average

# Print the new matrix
print(new_matrix)
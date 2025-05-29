# Define the original matrix
matrix = [
    [-8, 17, 7, 9],
    [12, -10, 10, -3],
    [12, 5, -6, 1],
    [-6, -10, 3, 15]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the new matrix
new_size = len(matrix) - window_size + 1

# Initialize the new matrix for storing averages
average_matrix = [[0] * new_size for _ in range(new_size)]

# Perform the average pooling operation
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
        average_matrix[i][j] = average

# Print the resulting average matrix
print(average_matrix)
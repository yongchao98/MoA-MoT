# Define the original matrix
matrix = [
    [-8, -6, 2, 18],
    [3, 13, 20, 9],
    [-8, -6, 0, -4],
    [-10, -4, 10, 10]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the resulting matrix
result_size = len(matrix) - window_size + 1

# Initialize the result matrix
result_matrix = [[0] * result_size for _ in range(result_size)]

# Perform average pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the 2x2 window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Calculate the average
        average = sum(window) / len(window)
        # Store the average in the result matrix
        result_matrix[i][j] = average

# Print the result matrix
print(result_matrix)
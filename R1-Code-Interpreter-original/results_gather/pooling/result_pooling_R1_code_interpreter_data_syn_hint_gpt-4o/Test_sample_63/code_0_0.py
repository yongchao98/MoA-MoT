# Define the original matrix
matrix = [
    [20, 8, 7, -1],
    [12, 3, -2, -3],
    [5, 4, 18, -10],
    [2, 10, -7, 6]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the new matrix
new_size = len(matrix) - window_size + 1

# Initialize the new matrix to store the averages
result_matrix = [[0] * new_size for _ in range(new_size)]

# Perform the average pooling
for i in range(new_size):
    for j in range(new_size):
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
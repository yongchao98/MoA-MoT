# Define the original matrix
matrix = [
    [11, 15, 18, 18],
    [20, -8, -4, 19],
    [14, -9, 11, -7],
    [-9, 2, -7, 2]
]

# Define the size of the sliding window
window_size = 2

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = [[0] * output_size for _ in range(output_size)]

# Perform max pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the current 2x2 window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Find the maximum value in the current window
        max_value = max(window)
        # Assign the maximum value to the output matrix
        output_matrix[i][j] = max_value

# Print the output matrix
print(output_matrix)
# Define the original matrix
matrix = [
    [5, 18, 17],
    [16, 7, 5],
    [-9, 4, 3]
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
        # Extract the current window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Find the maximum value in the current window
        output_matrix[i][j] = max(window)

# Print the output matrix
print(output_matrix)
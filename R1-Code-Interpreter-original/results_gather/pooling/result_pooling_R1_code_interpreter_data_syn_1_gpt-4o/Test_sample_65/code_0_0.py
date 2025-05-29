# Define the original matrix
matrix = [
    [16, 2, -9, 17],
    [5, -3, -5, 2],
    [3, 7, 19, -9],
    [8, 5, -8, 12]
]

# Define the size of the sliding window
window_size = 2

# Calculate the dimensions of the resulting matrix
result_rows = len(matrix) - window_size + 1
result_cols = len(matrix[0]) - window_size + 1

# Initialize the resulting matrix
result_matrix = [[0] * result_cols for _ in range(result_rows)]

# Perform the average pooling operation
for i in range(result_rows):
    for j in range(result_cols):
        # Extract the 2x2 window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Calculate the average of the window
        window_average = sum(window) / len(window)
        # Store the average in the result matrix
        result_matrix[i][j] = window_average

# Print the resulting matrix
print(result_matrix)
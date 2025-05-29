# Define the input matrix
matrix = [
    [17, -6, 2, 10, 6, 1],
    [5, -8, -4, 12, 7, -2],
    [10, 6, -2, -3, 4, 2],
    [6, 5, 13, -2, 18, 15],
    [10, -8, 6, 5, -9, -8],
    [17, -7, 16, 7, -4, 7]
]

# Define the size of the sliding window
window_size = 3

# Calculate the dimensions of the output matrix
output_rows = len(matrix) - window_size + 1
output_cols = len(matrix[0]) - window_size + 1

# Initialize the output matrix
output_matrix = [[0] * output_cols for _ in range(output_rows)]

# Perform average pooling
for i in range(output_rows):
    for j in range(output_cols):
        # Extract the current window
        window = [matrix[i+x][j+y] for x in range(window_size) for y in range(window_size)]
        # Calculate the average of the current window
        window_average = sum(window) / (window_size * window_size)
        # Store the average in the output matrix
        output_matrix[i][j] = window_average

# Print the output matrix
print(output_matrix)
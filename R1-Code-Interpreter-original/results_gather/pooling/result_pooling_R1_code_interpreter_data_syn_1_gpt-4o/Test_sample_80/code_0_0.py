# Define the input matrix
matrix = [
    [0, 5, 3, 10, -7],
    [-10, -3, 11, -3, 19],
    [-3, 2, 5, -5, -1],
    [2, -9, 12, 8, 5],
    [3, 18, -6, 6, 15]
]

# Define the size of the sliding window
window_size = 3

# Calculate the dimensions of the output matrix
output_rows = len(matrix) - window_size + 1
output_cols = len(matrix[0]) - window_size + 1

# Initialize the output matrix
output_matrix = []

# Perform the average pooling
for i in range(output_rows):
    row = []
    for j in range(output_cols):
        # Extract the current window
        window = [matrix[i+x][j+y] for x in range(window_size) for y in range(window_size)]
        # Calculate the average of the current window
        average = sum(window) / (window_size * window_size)
        row.append(average)
    output_matrix.append(row)

# Print the output matrix
print(output_matrix)
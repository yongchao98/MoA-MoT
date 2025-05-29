# Define the input matrix
matrix = [
    [-6, 19, 11, 5, -6, 20],
    [6, -9, 10, 9, 0, 7],
    [15, -7, 17, 18, 0, 16],
    [2, 1, 0, 6, 6, 10],
    [15, -10, 20, -4, -8, 16],
    [3, 11, -5, 6, 16, 12]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the result matrix
result = [[0] * output_size for _ in range(output_size)]

# Perform average pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the current window
        window = [matrix[i+x][j+y] for x in range(window_size) for y in range(window_size)]
        # Calculate the average of the window
        average = sum(window) / (window_size * window_size)
        # Store the average in the result matrix
        result[i][j] = average

# Print the result matrix
print(result)
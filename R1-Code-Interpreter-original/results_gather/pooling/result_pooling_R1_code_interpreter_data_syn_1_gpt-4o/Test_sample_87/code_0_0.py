# Define the original matrix
matrix = [
    [2, -6, -2, 17, 3, 1],
    [-4, 0, -10, 1, 3, -4],
    [20, -4, 3, 12, -9, -2],
    [-2, -7, -2, 19, -1, 1],
    [10, 11, 4, 0, 5, 17],
    [-4, 0, -7, 9, -8, 6]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the output matrix
output_size = len(matrix) - window_size + 1

# Initialize the output matrix
output_matrix = [[0] * output_size for _ in range(output_size)]

# Perform average pooling
for i in range(output_size):
    for j in range(output_size):
        # Extract the current 3x3 window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Flatten the window and calculate the average
        window_values = [num for row in window for num in row]
        average = sum(window_values) / len(window_values)
        # Store the average in the output matrix
        output_matrix[i][j] = average

# Print the output matrix
print(output_matrix)
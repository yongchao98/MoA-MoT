# Define the original 6x6 matrix
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

# Calculate the size of the resulting matrix
result_size = len(matrix) - window_size + 1

# Initialize the result matrix
result_matrix = [[0] * result_size for _ in range(result_size)]

# Perform average pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the 3x3 window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Flatten the window and calculate the average
        window_values = [num for row in window for num in row]
        average = sum(window_values) / len(window_values)
        # Store the average in the result matrix
        result_matrix[i][j] = average

# Print the result matrix
print(result_matrix)
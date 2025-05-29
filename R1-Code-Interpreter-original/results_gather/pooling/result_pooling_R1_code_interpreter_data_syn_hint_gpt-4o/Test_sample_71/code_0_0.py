# Define the original matrix
matrix = [
    [0, 12, 3, 2, -9],
    [-1, -1, -7, 10, 1],
    [20, -8, -3, -10, -6],
    [12, -4, 5, 3, 13],
    [12, 12, 15, 3, 0]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the resulting matrix after pooling
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
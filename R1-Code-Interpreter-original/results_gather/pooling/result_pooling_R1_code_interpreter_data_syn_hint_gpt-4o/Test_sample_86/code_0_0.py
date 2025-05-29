# Define the original matrix
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

# Calculate the size of the result matrix
result_size = len(matrix) - window_size + 1

# Initialize the result matrix
result_matrix = [[0] * result_size for _ in range(result_size)]

# Perform average pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the current 3x3 window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Flatten the window and calculate the average
        window_values = [num for row in window for num in row]
        average = sum(window_values) / len(window_values)
        # Store the average in the result matrix
        result_matrix[i][j] = average

# Print the result matrix
print(result_matrix)
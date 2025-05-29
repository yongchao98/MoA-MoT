# Define the 5x5 matrix
matrix = [
    [15, 5, 7, 1, -6],
    [-10, 19, 14, 0, 0],
    [-5, 16, 15, -9, 16],
    [12, -3, 3, 9, 12],
    [9, 7, 8, 16, -2]
]

# Size of the sliding window
window_size = 3

# Calculate the size of the result matrix
result_size = len(matrix) - window_size + 1

# Initialize the result matrix
result_matrix = [[0] * result_size for _ in range(result_size)]

# Perform average pooling
for i in range(result_size):
    for j in range(result_size):
        # Extract the 3x3 window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Flatten the window and calculate the average
        window_elements = [element for row in window for element in row]
        average = sum(window_elements) / len(window_elements)
        # Store the average in the result matrix
        result_matrix[i][j] = average

# Print the result matrix
print(result_matrix)
# Define the original matrix
matrix = [
    [20, -6, 8, 6, -5],
    [0, 4, 19, 2, 18],
    [4, 1, 19, 15, -3],
    [14, 5, 16, 19, -6],
    [3, 3, -1, -7, 20]
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
        # Extract the current window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Flatten the window and calculate the average
        window_values = [num for row in window for num in row]
        average = sum(window_values) / len(window_values)
        # Store the average in the output matrix
        output_matrix[i][j] = average

# Print the output matrix
print(output_matrix)
# Define the original matrix
matrix = [
    [8, 13, 12, 13, -5],
    [16, 7, 6, 12, 15],
    [-4, -3, -5, 13, 1],
    [0, 5, 12, 6, -2],
    [-5, 17, -10, -8, 7]
]

# Define the size of the sliding window
window_size = 3

# Calculate the size of the new matrix
new_size = len(matrix) - window_size + 1

# Initialize the new matrix
new_matrix = [[0] * new_size for _ in range(new_size)]

# Perform max pooling
for i in range(new_size):
    for j in range(new_size):
        # Extract the current window
        window = [matrix[x][j:j+window_size] for x in range(i, i+window_size)]
        # Find the maximum value in the current window
        max_value = max(max(row) for row in window)
        # Assign the maximum value to the new matrix
        new_matrix[i][j] = max_value

# Print the new matrix
print(new_matrix)
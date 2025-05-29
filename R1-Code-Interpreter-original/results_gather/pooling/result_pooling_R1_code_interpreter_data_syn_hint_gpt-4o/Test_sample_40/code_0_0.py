# Original matrix
matrix = [
    [-4, -9, -1, 7, -8],
    [13, 0, 8, 8, 9],
    [0, 19, -6, 0, 19],
    [7, 12, 12, -1, -3],
    [9, -5, 5, 11, 0]
]

# Dimensions
N = 5
n = 3

# Resultant matrix after max pooling
result = []

# Perform max pooling
for i in range(N - n + 1):
    row = []
    for j in range(N - n + 1):
        # Extract the 3x3 window
        window = [matrix[x][j:j+n] for x in range(i, i+n)]
        # Find the maximum value in the window
        max_value = max(max(row) for row in window)
        row.append(max_value)
    result.append(row)

# Print the result matrix
print(result)
# Original matrix
matrix = [
    [-7, 19, -9, -8, 6],
    [6, 11, 4, -6, -10],
    [0, 1, 6, 16, -9],
    [-10, 13, 14, -8, -1],
    [12, 18, 10, 20, 5]
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

# Print the result
print(result)
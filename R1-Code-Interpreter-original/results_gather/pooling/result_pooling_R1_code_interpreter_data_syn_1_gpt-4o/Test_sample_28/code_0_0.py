# Original matrix
matrix = [
    [9, -9, 4, -3],
    [-1, -8, 15, 7],
    [11, 13, -9, 17],
    [20, 15, 20, -7]
]

# Dimensions of the original matrix
N = 4

# Size of the sliding window
n = 2

# Resultant matrix after max pooling
result = []

# Perform max pooling
for i in range(N - n + 1):
    row = []
    for j in range(N - n + 1):
        # Extract the 2x2 window
        window = [
            matrix[i][j], matrix[i][j+1],
            matrix[i+1][j], matrix[i+1][j+1]
        ]
        # Find the maximum value in the window
        max_value = max(window)
        row.append(max_value)
    result.append(row)

# Print the result matrix
print(result)
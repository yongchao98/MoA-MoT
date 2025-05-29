# Original matrix
matrix = [
    [6, -6, -9, -8],
    [2, 7, 4, 16],
    [5, 19, -6, 11],
    [20, 14, -8, 4]
]

# Dimensions of the original matrix
N = 4

# Dimensions of the sliding window
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

# Print the result
print(result)
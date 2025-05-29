# Original matrix
matrix = [
    [4, -2, -2, -3],
    [-8, 17, 1, 19],
    [15, 14, -9, 7],
    [-9, 3, 3, 19]
]

# Dimensions of the original matrix
N = 4
# Dimensions of the sliding window
n = 2

# Resultant matrix after max pooling
result = []

# Iterate over the matrix with a sliding window
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
        # Append the max value to the current row
        row.append(max_value)
    # Append the row to the result matrix
    result.append(row)

# Print the result matrix
print(result)
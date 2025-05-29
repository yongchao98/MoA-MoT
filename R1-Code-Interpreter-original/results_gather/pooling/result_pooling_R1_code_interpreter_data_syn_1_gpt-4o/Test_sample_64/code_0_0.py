# Original matrix
matrix = [
    [7, 18, 18, -4],
    [12, 6, 14, -6],
    [1, -2, -8, 3],
    [20, 16, 9, -6]
]

# Dimensions of the original matrix
N = 4

# Size of the sliding window
n = 2

# Resultant matrix after average pooling
result = []

# Perform average pooling
for i in range(N - n + 1):
    row = []
    for j in range(N - n + 1):
        # Calculate the average of the current 2x2 window
        window_sum = (matrix[i][j] + matrix[i][j+1] +
                      matrix[i+1][j] + matrix[i+1][j+1])
        window_avg = window_sum / (n * n)
        row.append(window_avg)
    result.append(row)

# Print the result
print(result)
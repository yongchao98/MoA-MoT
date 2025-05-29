# Original matrix
matrix = [
    [-8, -6, 2, 18],
    [3, 13, 20, 9],
    [-8, -6, 0, -4],
    [-10, -4, 10, 10]
]

# Dimensions of the original matrix
N = 4

# Dimensions of the sliding window
n = 2

# Resultant matrix after average pooling
result = []

# Iterate over the matrix with the sliding window
for i in range(N - n + 1):
    row = []
    for j in range(N - n + 1):
        # Calculate the average of the current 2x2 window
        window_sum = (matrix[i][j] + matrix[i][j+1] +
                      matrix[i+1][j] + matrix[i+1][j+1])
        window_avg = window_sum / (n * n)
        row.append(window_avg)
    result.append(row)

# Print the resultant matrix
print(result)
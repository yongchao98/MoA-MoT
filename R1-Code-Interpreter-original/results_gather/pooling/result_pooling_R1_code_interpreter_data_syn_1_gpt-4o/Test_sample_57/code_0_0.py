# Define the original matrix
matrix = [
    [17, 8, 5, 8],
    [10, -6, -6, 4],
    [0, 4, 16, -9],
    [17, 7, 6, -8]
]

# Dimensions of the original matrix
N = 4

# Size of the sliding window
n = 2

# Initialize the result matrix for average pooling
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

# Print the result matrix
print(result)
# Define the original matrix
matrix = [
    [16, 2, -9, 17],
    [5, -3, -5, 2],
    [3, 7, 19, -9],
    [8, 5, -8, 12]
]

# Dimensions of the original matrix
N = 4
# Dimensions of the sliding window
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
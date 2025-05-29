# Define the original matrix
matrix = [
    [10, 8, -7, -3, 19],
    [-8, -5, 14, 11, -2],
    [14, -5, 9, 8, -9],
    [2, 10, -5, 7, 12],
    [18, 0, -9, -10, 1]
]

# Dimensions of the matrix
N = 5
n = 3

# Resultant matrix after average pooling
result = []

# Iterate over the matrix with a sliding window of size 3x3
for i in range(N - n + 1):
    row = []
    for j in range(N - n + 1):
        # Calculate the sum of the current 3x3 window
        window_sum = 0
        for k in range(n):
            for l in range(n):
                window_sum += matrix[i + k][j + l]
        # Calculate the average
        window_avg = window_sum / (n * n)
        row.append(window_avg)
    result.append(row)

# Print the result
print(result)
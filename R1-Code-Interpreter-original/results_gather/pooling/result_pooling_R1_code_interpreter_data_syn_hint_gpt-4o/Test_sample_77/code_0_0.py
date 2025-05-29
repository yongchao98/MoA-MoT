# Original 5x5 matrix
matrix = [
    [-3, -4, -3, -5, 5],
    [-6, 13, 10, 13, 5],
    [20, 15, -6, 8, 9],
    [8, 9, 4, 4, 10],
    [-1, 7, -2, 1, -7]
]

# Dimensions
N = 5
n = 3

# Resultant matrix for average pooling
result = []

# Perform average pooling
for i in range(N - n + 1):
    row = []
    for j in range(N - n + 1):
        # Calculate the sum of the 3x3 window
        window_sum = 0
        for x in range(n):
            for y in range(n):
                window_sum += matrix[i + x][j + y]
        # Calculate the average
        window_avg = window_sum / (n * n)
        row.append(window_avg)
    result.append(row)

# Print the result
print(result)
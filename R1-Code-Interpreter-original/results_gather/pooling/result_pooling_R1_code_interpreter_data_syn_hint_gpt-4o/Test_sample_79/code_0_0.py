# Original matrix
matrix = [
    [20, -6, 8, 6, -5],
    [0, 4, 19, 2, 18],
    [4, 1, 19, 15, -3],
    [14, 5, 16, 19, -6],
    [3, 3, -1, -7, 20]
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
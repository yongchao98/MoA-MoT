# Original 5x5 matrix
matrix = [
    [15, 5, 7, 1, -6],
    [-10, 19, 14, 0, 0],
    [-5, 16, 15, -9, 16],
    [12, -3, 3, 9, 12],
    [9, 7, 8, 16, -2]
]

# Dimensions
N = 5
n = 3

# Resultant matrix after average pooling
result = []

# Perform average pooling
for i in range(N - n + 1):
    row = []
    for j in range(N - n + 1):
        # Calculate the sum of the current 3x3 window
        window_sum = 0
        for x in range(n):
            for y in range(n):
                window_sum += matrix[i + x][j + y]
        # Calculate the average
        window_average = window_sum / (n * n)
        row.append(window_average)
    result.append(row)

# Print the result
print(result)
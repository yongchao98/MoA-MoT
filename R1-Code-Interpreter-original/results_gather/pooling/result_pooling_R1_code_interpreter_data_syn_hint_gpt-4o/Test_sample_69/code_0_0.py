# Define the original matrix
matrix = [
    [8, 11, 19, 19, 0],
    [7, 11, 15, 8, 10],
    [11, 7, 4, 15, 7],
    [4, -6, -9, 10, 12],
    [-3, 15, 20, 9, 15]
]

# Dimensions of the matrix and the sliding window
N = 5
n = 3

# Initialize the result matrix
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
        window_avg = window_sum / (n * n)
        row.append(window_avg)
    result.append(row)

# Print the result matrix
print(result)
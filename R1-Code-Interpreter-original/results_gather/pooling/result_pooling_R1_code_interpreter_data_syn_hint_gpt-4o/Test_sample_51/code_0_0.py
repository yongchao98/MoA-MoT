# Define the original matrix
matrix = [
    [10, 11, -2, 1],
    [5, 18, 18, 18],
    [-9, -1, 7, 11],
    [12, -1, -4, 14]
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
        # Calculate the sum of the current 2x2 window
        window_sum = (matrix[i][j] + matrix[i][j+1] +
                      matrix[i+1][j] + matrix[i+1][j+1])
        # Calculate the average
        window_avg = window_sum / (n * n)
        # Append the average to the current row
        row.append(window_avg)
    # Append the row to the result matrix
    result.append(row)

# Print the result matrix
print(result)
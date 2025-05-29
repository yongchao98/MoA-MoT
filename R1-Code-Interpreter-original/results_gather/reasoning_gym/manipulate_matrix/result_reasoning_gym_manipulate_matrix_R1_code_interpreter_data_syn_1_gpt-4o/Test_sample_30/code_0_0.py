# Define the original matrix
matrix = [
    [3, 7, 4],
    [5, 2, 8],
    [2, 2, 0],
    [4, 2, 1]
]

# Get the dimensions of the matrix
n = len(matrix)    # Number of rows
m = len(matrix[0]) # Number of columns

# Create a new matrix to store the result of the counterdiagonal mirroring
mirrored_matrix = [[0] * n for _ in range(m)]

# Perform the counterdiagonal mirroring
for i in range(n):
    for j in range(m):
        mirrored_matrix[m-j-1][n-i-1] = matrix[i][j]

# Print the final mirrored matrix
for row in mirrored_matrix:
    print(' '.join(map(str, row)))
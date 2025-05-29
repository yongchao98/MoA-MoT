# Define the input matrix as a list of lists
matrix = [
    [8, 9, 0, 4, 2, 3, 4],
    [1, 3, 8, 0, 9, 3, 8],
    [0, 8, 5, 7, 2, 5, 2],
    [0, 5, 9, 0, 6, 9, 3],
    [2, 5, 6, 7, 1, 0, 3],
    [9, 2, 2, 3, 8, 5, 4],
    [9, 0, 4, 5, 3, 6, 0]
]

# Get the size of the matrix
n = len(matrix)

# Create a new matrix for the rotated result
rotated_matrix = [[0] * n for _ in range(n)]

# Perform the rotation
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n - 1 - i] = matrix[i][j]

# Format the output
output = "\n".join(" ".join(map(str, row)) for row in rotated_matrix)
print(f"<<<{output}>>>")
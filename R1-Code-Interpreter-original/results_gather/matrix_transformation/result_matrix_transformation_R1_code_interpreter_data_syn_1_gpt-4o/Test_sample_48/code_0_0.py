# Original matrix
matrix = [
    ["F", "M", "R"],
    ["M", "Q", "W"],
    ["H", "S", "C"]
]

# Dimensions of the matrix
n = len(matrix)

# Create a new matrix for the rotated result
rotated_matrix = [['' for _ in range(n)] for _ in range(n)]

# Rotate the matrix 90 degrees clockwise
for i in range(n):
    for j in range(n):
        rotated_matrix[j][n-1-i] = matrix[i][j]

# Flatten the rotated matrix and format the output
flattened_result = [item for row in rotated_matrix for item in row]
print(flattened_result)
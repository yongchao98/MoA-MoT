# Original matrix
matrix = [
    ["8", "3", "2"],
    ["6", "3", "1"]
]

# Number of rows and columns
n = len(matrix)
m = len(matrix[0])

# Create a new matrix for the rotated result
rotated_matrix = [["" for _ in range(n)] for _ in range(m)]

# Perform the rotation
for i in range(n):
    for j in range(m):
        rotated_matrix[j][n-1-i] = matrix[i][j]

# Flatten the rotated matrix and format the output
flattened_result = [item for sublist in rotated_matrix for item in sublist]
print(flattened_result)
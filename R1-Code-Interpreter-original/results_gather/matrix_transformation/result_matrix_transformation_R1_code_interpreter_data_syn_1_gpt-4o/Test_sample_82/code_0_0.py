# Define the original matrix
matrix = [
    ['*', '&', '&', '@'],
    ['&', '%', '*', '!'],
    ['!', '*', '%', '^'],
    ['^', '@', '%', '$']
]

# Get the number of rows (n) and columns (m)
n = len(matrix)
m = len(matrix[0])

# Initialize the rotated matrix with dimensions m x n
rotated_matrix = [['' for _ in range(n)] for _ in range(m)]

# Perform the rotation
for i in range(n):
    for j in range(m):
        rotated_matrix[j][n - 1 - i] = matrix[i][j]

# Flatten the rotated matrix to a single list
flattened_result = [item for row in rotated_matrix for item in row]

# Print the result
print(flattened_result)
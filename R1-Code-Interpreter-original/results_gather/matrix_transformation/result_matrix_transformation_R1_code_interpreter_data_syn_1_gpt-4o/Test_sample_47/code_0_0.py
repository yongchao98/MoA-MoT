# Original matrix
matrix = [
    ['S', 'P', 'W'],
    ['O', 'D', 'A'],
    ['M', 'J', 'Y']
]

# Number of rows and columns
n = len(matrix)
m = len(matrix[0])

# Initialize the rotated matrix
rotated = [['' for _ in range(n)] for _ in range(m)]

# Fill the rotated matrix
for i in range(n):
    for j in range(m):
        rotated[j][n-1-i] = matrix[i][j]

# Flatten the rotated matrix to a single list
flattened_rotated = [item for sublist in rotated for item in sublist]

# Print the result in the required format
print(flattened_rotated)
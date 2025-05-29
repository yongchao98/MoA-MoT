# Original matrix
matrix = [
    ['6', 'P', 'T', 'D'],
    ['7', 'L', 'L', '4'],
    ['P', 'C', 'O', 'S']
]

# Number of rows and columns
n = len(matrix)
m = len(matrix[0])

# Initialize the rotated matrix
rotated_matrix = [['' for _ in range(n)] for _ in range(m)]

# Fill the rotated matrix
for i in range(n):
    for j in range(m):
        rotated_matrix[j][n-1-i] = matrix[i][j]

# Flatten the rotated matrix to a list of values
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)
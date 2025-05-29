# Original matrix
matrix = [
    ['L', '1', '9', 'M'],
    ['T', 'J', '9', 'G'],
    ['Z', 'X', 'U', 'U']
]

# Number of rows and columns
n = len(matrix)
m = len(matrix[0])

# Rotated matrix initialization
rotated_matrix = [['' for _ in range(n)] for _ in range(m)]

# Rotate the matrix 90 degrees clockwise
for i in range(n):
    for j in range(m):
        rotated_matrix[j][n-1-i] = matrix[i][j]

# Flatten the rotated matrix to a single list
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(f"<<<{','.join(flattened_result)}>>>")
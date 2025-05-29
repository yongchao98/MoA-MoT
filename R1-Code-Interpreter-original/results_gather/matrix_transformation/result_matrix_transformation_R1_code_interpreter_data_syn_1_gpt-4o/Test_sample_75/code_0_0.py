# Original matrix
matrix = [
    ['N', '8', 'C', '3'],
    ['C', 'L', '2', 'L'],
    ['N', 'B', 'L', 'W']
]

# Number of rows and columns
rows = len(matrix)
cols = len(matrix[0])

# Initialize the rotated matrix
rotated_matrix = [['' for _ in range(rows)] for _ in range(cols)]

# Fill the rotated matrix
for i in range(cols):
    for j in range(rows):
        rotated_matrix[i][j] = matrix[rows - j - 1][i]

# Flatten the rotated matrix for output
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)
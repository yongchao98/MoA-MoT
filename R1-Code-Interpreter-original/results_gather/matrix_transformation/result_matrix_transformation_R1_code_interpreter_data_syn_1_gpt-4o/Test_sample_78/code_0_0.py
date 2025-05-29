# Original matrix
matrix = [
    ['L', '1', '9', 'M'],
    ['T', 'J', '9', 'G'],
    ['Z', 'X', 'U', 'U']
]

# Dimensions of the original matrix
rows = len(matrix)
cols = len(matrix[0])

# Create a new matrix for the rotated result
rotated_matrix = [['' for _ in range(rows)] for _ in range(cols)]

# Fill the rotated matrix
for i in range(rows):
    for j in range(cols):
        rotated_matrix[j][rows - 1 - i] = matrix[i][j]

# Flatten the rotated matrix and format the output
flattened_result = [item for sublist in rotated_matrix for item in sublist]
print(flattened_result)
# Original matrix
matrix = [
    ['5', '9', '4'],
    ['8', '2', '1']
]

# Number of rows and columns in the original matrix
rows = len(matrix)
cols = len(matrix[0])

# Create a new matrix for the rotated result
rotated_matrix = [['' for _ in range(rows)] for _ in range(cols)]

# Apply the transformation
for i in range(rows):
    for j in range(cols):
        rotated_matrix[j][rows - 1 - i] = matrix[i][j]

# Flatten the rotated matrix to a list of values
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(f"<<<{','.join(flattened_result)}>>>")
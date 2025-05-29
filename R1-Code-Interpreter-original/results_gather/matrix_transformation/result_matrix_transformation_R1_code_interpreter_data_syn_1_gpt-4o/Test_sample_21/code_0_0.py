# Original matrix
matrix = [
    [2, 1, 3],
    [3, 7, 7]
]

# Number of rows and columns in the original matrix
rows = len(matrix)
cols = len(matrix[0])

# Initialize the rotated matrix with dimensions cols x rows
rotated_matrix = [[0] * rows for _ in range(cols)]

# Fill the rotated matrix
for i in range(rows):
    for j in range(cols):
        rotated_matrix[j][rows - 1 - i] = matrix[i][j]

# Flatten the rotated matrix to a list of values
flattened_result = [str(item) for sublist in rotated_matrix for item in sublist]

# Print the result in the required format
print(f"<<<{','.join(flattened_result)}>>>")
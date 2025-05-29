# Original matrix
matrix = [
    ["M", "C", "K", "X"],
    ["8", "7", "N", "C"],
    ["Q", "M", "V", "1"]
]

# Number of rows and columns in the original matrix
rows = len(matrix)
cols = len(matrix[0])

# Initialize the rotated matrix with dimensions cols x rows
rotated_matrix = [['' for _ in range(rows)] for _ in range(cols)]

# Fill the rotated matrix
for i in range(rows):
    for j in range(cols):
        rotated_matrix[j][rows - 1 - i] = matrix[i][j]

# Flatten the rotated matrix for output
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)
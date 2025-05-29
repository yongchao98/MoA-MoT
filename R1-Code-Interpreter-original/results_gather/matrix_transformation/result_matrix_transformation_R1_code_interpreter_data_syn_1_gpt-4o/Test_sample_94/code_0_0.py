# Original matrix
matrix = [
    ['*', '#', '^'],
    ['!', '*', '*', '&'],
    ['#', '*', '%', '@'],
    ['%', '@', '%', '#']
]

# Determine the number of rows and columns in the original matrix
num_rows = len(matrix)
num_cols = len(matrix[0])

# Create a new matrix with swapped dimensions
rotated_matrix = [['' for _ in range(num_rows)] for _ in range(num_cols)]

# Populate the new matrix with rotated values
for i in range(num_rows):
    for j in range(num_cols):
        rotated_matrix[j][num_rows - 1 - i] = matrix[i][j]

# Flatten the rotated matrix for output
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)
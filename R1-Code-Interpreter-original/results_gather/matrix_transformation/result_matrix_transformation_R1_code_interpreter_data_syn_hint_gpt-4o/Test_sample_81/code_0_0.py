# Original matrix
matrix = [
    ['&', '^', '!', '$'],
    ['#', '&', '%', '!'],
    ['#', '#', '&', '&'],
    ['#', '%', '$', '!']
]

# Number of rows and columns
n = len(matrix)
m = len(matrix[0])

# Rotated matrix initialization
rotated_matrix = []

# Rotate the matrix 90 degrees clockwise
for col in range(m):
    new_row = []
    for row in range(n-1, -1, -1):
        new_row.append(matrix[row][col])
    rotated_matrix.append(new_row)

# Flatten the rotated matrix for output
flattened_result = [item for sublist in rotated_matrix for item in sublist]

# Print the result
print(flattened_result)
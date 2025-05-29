# Original matrix
matrix = [
    ['R', '7', 'I', 'S'],
    ['S', 'K', 'C', 'A'],
    ['1', '1', 'S', 'E']
]

# Number of rows and columns
rows = len(matrix)
cols = len(matrix[0])

# Initialize the rotated matrix
rotated_matrix = []

# Rotate the matrix 90 degrees clockwise
for col in range(cols):
    new_row = []
    for row in range(rows - 1, -1, -1):
        new_row.append(matrix[row][col])
    rotated_matrix.append(new_row)

# Flatten the rotated matrix and format as a single string
flattened_result = [item for sublist in rotated_matrix for item in sublist]
formatted_result = ','.join(flattened_result)

# Print the result
print(formatted_result)
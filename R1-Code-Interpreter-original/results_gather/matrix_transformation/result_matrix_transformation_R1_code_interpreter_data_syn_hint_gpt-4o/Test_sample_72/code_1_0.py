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

# Format the rotated matrix as a list of strings
formatted_result = [','.join(row) for row in rotated_matrix]

# Print the result
print(formatted_result)
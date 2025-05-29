# Original matrix
matrix = [
    ['G', '6', 'B', 'H'],
    ['K', 'Q', '3', '7'],
    ['7', '2', 'H', '4']
]

# Transpose the matrix
transposed_matrix = list(zip(*matrix))

# Reverse each row of the transposed matrix
rotated_matrix = [list(row)[::-1] for row in transposed_matrix]

# Print the rotated matrix
print(rotated_matrix)
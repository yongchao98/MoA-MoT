# Original matrix
matrix = [
    [6, 2, 8],
    [2, 9, 7]
]

# Transpose the matrix
transposed_matrix = list(zip(*matrix))

# Reverse each row of the transposed matrix
rotated_matrix = [list(row)[::-1] for row in transposed_matrix]

# Print the rotated matrix
print(rotated_matrix)
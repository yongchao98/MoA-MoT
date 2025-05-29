def rotate_90_clockwise(matrix):
    # Transpose the matrix
    transposed_matrix = list(zip(*matrix))
    # Reverse each row to get the rotated matrix
    rotated_matrix = [list(reversed(row)) for row in transposed_matrix]
    return rotated_matrix

# Original matrix
matrix = [['R', 'L'], ['K', 'N']]

# Rotate the matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix in the required format
print(rotated_matrix)
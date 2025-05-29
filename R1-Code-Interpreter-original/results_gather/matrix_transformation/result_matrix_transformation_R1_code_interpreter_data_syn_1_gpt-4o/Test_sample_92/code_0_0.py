# Original matrix
matrix = [
    ['*', '#', '$', '*'],
    ['%', '@', '#', '@'],
    ['&', '*', '^', '%'],
    ['^', '!', '$', '*']
]

# Rotate the matrix 90 degrees clockwise
def rotate_90_clockwise(mat):
    # Transpose the matrix
    transposed = list(zip(*mat))
    # Reverse each row of the transposed matrix
    rotated = [list(row)[::-1] for row in transposed]
    return rotated

# Get the rotated matrix
rotated_matrix = rotate_90_clockwise(matrix)

# Print the rotated matrix
print(rotated_matrix)
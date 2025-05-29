def rotate_180(matrix):
    # Reverse the order of rows and then reverse each row
    return [row[::-1] for row in matrix[::-1]]

# Input matrix
matrix = [
    [2, 2, 3, 6, 2],
    [8, 6, 1, 8, 9],
    [5, 5, 7, 0, 2],
    [5, 6, 1, 2, 6],
    [6, 6, 6, 0, 2]
]

# Rotate the matrix by 180 degrees
rotated_matrix = rotate_180(matrix)

# Print the rotated matrix
print(rotated_matrix)
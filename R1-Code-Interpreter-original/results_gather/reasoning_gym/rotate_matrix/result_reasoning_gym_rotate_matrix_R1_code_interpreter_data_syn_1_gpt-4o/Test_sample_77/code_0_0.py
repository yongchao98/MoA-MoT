def rotate_180(matrix):
    # Reverse the order of rows and then reverse each row
    return [row[::-1] for row in matrix[::-1]]

# Input matrix
matrix = [
    [3, 7, 7, 8, 0, 1],
    [3, 8, 7, 4, 9, 8],
    [1, 3, 7, 8, 4, 6],
    [2, 4, 8, 8, 5, 1],
    [8, 2, 0, 5, 9, 3],
    [2, 9, 2, 6, 7, 2]
]

# Rotate the matrix by 180 degrees
rotated_matrix = rotate_180(matrix)

# Print the rotated matrix
print(rotated_matrix)
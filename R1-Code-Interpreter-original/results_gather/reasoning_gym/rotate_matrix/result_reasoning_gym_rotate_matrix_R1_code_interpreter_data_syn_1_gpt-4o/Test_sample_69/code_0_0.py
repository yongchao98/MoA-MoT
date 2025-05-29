def rotate_180(matrix):
    # Reverse the order of rows and then reverse each row
    return [row[::-1] for row in matrix[::-1]]

# Given matrix
matrix = [
    [8, 5, 7, 0, 8],
    [3, 9, 3, 4, 7],
    [3, 5, 7, 1, 4],
    [3, 2, 1, 7, 8],
    [5, 3, 7, 5, 3]
]

# Rotate the matrix by 180 degrees
rotated_matrix = rotate_180(matrix)

# Print the rotated matrix
print(rotated_matrix)